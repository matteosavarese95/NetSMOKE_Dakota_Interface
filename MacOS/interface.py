import sys
import os
from genericpath import isfile
import numpy as np
import pandas as pd

# ---------------- USER DEFINED PARAMETERS ---------------- #
P = 48.0        # kW
y_nh3 = 1.0     # mole fraction

# ---------------- HELPER FUNCTIONS ---------------- #
# calculate fuel and air mass flowrates
def calc_inlet_mass(Power, y_nh3, phi):
    
    # This function calculates the mass flow rate of fuel
    # and total air (based on the global equivalence ratio)
    # given the power (kW), the NH3 mole fraction of the fuel and 
    # The global equivalence ratio.

    # Calculate the stoichiometric coefficients
    y_h2 = 1 - y_nh3
    y = y_h2*2 + y_nh3*3
    w = y_nh3
    x = 0.
    m = x + y/4                                 # O2/fuel molar basis
    f = 0.79/0.21                               # Ratio between N2 and O2 in air
    af_mol  = m + m*f                           # Air/fuel stoichiometric molar basis
    af_mass = af_mol*29/(y_h2*2 + y_nh3*17)     # Air/fuel stoichiometric mass basis
    air_ecc = 1/phi -1                          # Air excess
    af_eff  = af_mass*(1+air_ecc)               # Air/Fuel operative mass basis

    # Calculate the fuel flowrate
    nH2O_s = y/2;                               # Stoichiometric moles of H2O produced
    DH_r   = -241.86*nH2O_s + 46.01*y_nh3;      # DH of global reaction kJ/mol
    LHV = -1000*DH_r/(y_h2*2 + y_nh3*17);       # LHV kJ/kg
    m_fuel = Power/LHV;                         # Fuel flowrate kg/s

    # Calculate air flowrate
    m_air = m_fuel * af_eff                 # Air flowrate kg/s

    return m_fuel, m_air

# Calculate the mass flowrates of the network
def calc_mass_flowrates(Power, y_nh3, phi_rich, phi_lean):

    # This function calculates the mass flow rates across the reactors 
    # and stores them in a Matrix M[i,j] where M[i,j] is the mass flow rate
    # from reactor i to reactor j. Please note that this function defines
    # the internal structure of the reactor network. The user should modify
    # this function to change the reactor network model.
    # 
    m_fuel, m_air_tot = calc_inlet_mass(Power, y_nh3, phi_lean)
    m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)

    # Initialize the matrix of mass flowrates
    M = np.zeros((5,5))

    # Reactor 0 is a "Fake" reactor that represents the fuel inlet
    M[0,2] = m_fuel

    # Reactor 1 is a "Fake" reactor that represents the air (rich) inlet
    M[1,2] = m_air_rich

    # Reactor 2 is flame PSR where rich combustion occurs
    M[2,3] = m_fuel + m_air_rich

    # Reactor 3 is the fake inlet of the remaining air
    M[3,4] = m_air_tot

    # Reactor 4 is the final PFR
    M[2,4] = m_fuel + m_air_rich

    return M

# ---------------- READING DAKOTA INPUTS ---------------- #

# Get the name of the file which is passed as an argument
filename = sys.argv[1]

# Read the input file
with open(filename, 'r') as f:
    params = f.readlines()

# Get the number of parameters
n_params = int(params[0].split()[0])

# Extract their names and values and put them in two separate lists
names = []
values = []
counter = 1
while counter < n_params+1:
    names.append(params[counter].split()[1])
    values.append(params[counter].split()[0])
    counter += 1

# Extract the name of the outputs
n_output = int(params[counter].split()[0])
print('File has', n_output, 'outputs')

# Extract the names of the outputs
counter = counter + 1
pos = counter 
output_names = []
while counter < pos + n_output:
    s = params[counter].split()
    output_names.append(s[1].split(':')[1])
    counter += 1
print('There are', len(output_names), 'outputs')

# ---------------- SUBSTITUTING PARAMETERS INTO REACTORS ---------------- #
# Check if phi_rich or phi_lean are in the input name list
print(names)
if 'phi_rich' in names and 'phi_lean' in names:
    phi_rich = float(values[names.index('phi_rich')])
    phi_lean = float(values[names.index('phi_lean')])
else:
    print('Error: phi_rich or phi_lean not found in input file')
    sys.exit(1)

# Get mass flowrates
M = calc_mass_flowrates(P, y_nh3, phi_rich, phi_lean)

# Get the number of reactors and the list of files to open to replace strings
print('Checking which files to open...\n')
nr = 0
fname_list = []
for i in range(1000):
    if isfile('input.cstr.'+str(i)+'.dic'):
        fname_list.append('input.cstr.' + str(i) + '.dic')
        nr += 1
    elif isfile('input.pfr.'+str(i)+'.dic'):
        fname_list.append('input.pfr.' + str(i) + '.dic')
        nr += 1
    else:
        print('Found', nr, 'reactors')
        break

# Check if main input is present
if not isfile('input.dic'):
    print('Error: input.dic not found')
    sys.exit(1)

# Open and modify each file
for fname in fname_list:
    newlines = []    # New lines to write in the modified input
    counter = 0      # Counter to keep track of the line number
    with open(fname, 'r') as fid:
        lines = fid.readlines()
        for item in lines:
            line = item.split()
            # For each of the names check if it is in line
            for name in names:
                if name in line:
                    # If it is, replace the name with the value
                    index = line.index(name)
                    line[index] = values[names.index(name)]
            
            for subitem in line:
                if subitem =='M0':
                    line[line.index(subitem)] = str(M[0,2])
                elif subitem =='M1':
                    line[line.index(subitem)] = str(M[1,2])
                elif subitem =='M2':
                    line[line.index(subitem)] = str(M[2,4])
                elif subitem =='M3':
                    line[line.index(subitem)] = str(M[3,4])
                elif subitem =='M4':
                    line[line.index(subitem)] = str(M[2,4]+M[3,4])
        
            newlines.append(line)
            counter += 1

    # Write modified lines to file
    with open(fname, 'w') as f:
        for item in newlines:
            for ss in item:
                if ss == item[-1]:
                    f.write(ss + '\n')
                else:
                    f.write(ss + ' ')

print('New inputs are ready')
print('Now modifying the internal connections of the network...')

# ---------------- WRITING THE MAIN INPUT.DIC ---------------- #

# Modify input dictionary
with open('input.dic', 'r') as fid:
    lines = fid.readlines()
    newlines = []
    for item in lines:
        line = item.split()
        try:
            r1 = int(line[0])
            r2 = int(line[1])
            line[2] = str(M[r1,r2])
        except:
            pass

        for subitem in line:
            if subitem == 'M0':
                line[line.index(subitem)] = str(M[0,2])
            elif subitem == 'M1':
                line[line.index(subitem)] = str(M[1,2])
            elif subitem == 'M3':
                line[line.index(subitem)] = str(M[3,4])
            elif subitem == 'M4':
                line[line.index(subitem)] = str(M[2,4]+M[3,4])

    
        newlines.append(line)

with open('input.dic', 'w') as f:
    for item in newlines:
        for ss in item:
            if ss == item[-1]:
                f.write(ss + '\n')
            else:
                f.write(ss + ' ')

print('Main input file is ready')

# ---------------- RUNNING THE SIMULATION ---------------- #
os.system('sh Run.sh')

# ---------------- WRITING OUTPUTS ---------------- #
print('Reading outputs...')

output_values = []
for item in output_names:
    s = item.split('_')
    outname = 'Output/Reactor.' + s[2] + '/Output.out'
    dt = pd.read_csv(outname, sep = '\s+')

    H2O_out = dt['H2O_x(17)'].values[0]
    O2_out = dt['O2_x(15)'].values[0]
    
    if s[0] == 'T':
        output_values.append(dt['T[K](5)'].values[-1])
    elif s[0] == 'NO':
        output_values.append(dt['NO_x(21)'].values[-1]*1e6)
    elif s[0] == 'NH3':
        output_values.append(dt['NH3_x(32)'].values[-1]*1e6)

if len(output_values) != len(output_names):
    print('Error: output_values and output_names have different lengths')
    sys.exit(1)

# Write outputs to file
output_dakota = sys.argv[2]
with open(output_dakota, 'w') as f:
    for i in range(len(output_values)):
        f.write(str(output_values[i]) + ' ' + output_names[i] + '\n')


