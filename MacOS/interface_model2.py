import sys
import os
from genericpath import isfile
import numpy as np
import pandas as pd
import cantera as ct

# ---------------- USER DEFINED PARAMETERS ---------------- #
P = 48.0                # kW
y_nh3 = 1.0             # mole fraction
Tout  = 1200.0          # K
Power = 48.0            # kW
T_in_rich = 300.0       # K (fuel temperature)

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
    m_air_bypass = m_air_tot - m_air_rich

    M = np.zeros((8,8))

    # First combustion zone
    M[0,2] = m_fuel                 # To flame PSR
    M[1,2] = m_air_rich             # To flame PSR

    M[2,7] = m_fuel + m_air_rich    # From ignition to quenching PSR
    M[7,3] = m_fuel + m_air_rich    # From quenching to flame PFR
    M[3,4] = m_fuel + m_air_rich    # From flame to first stagnation reactor

    # Stagnation zone
    stagn1 = 0.1
    stagn2 = 0.15
    stagn3 = 1 - stagn1 - stagn2
    M[1,4] = m_air_bypass*stagn1
    M[1,5] = m_air_bypass*stagn2
    M[1,6] = m_air_bypass*stagn3

    M[4,5] = m_fuel + m_air_rich + m_air_bypass*stagn1
    M[5,6] = m_fuel + m_air_rich + m_air_bypass*stagn1 + m_air_bypass*stagn2


    return M

# ---------------- HELPER FUNCTIONS ---------------- #
# Calculate fuel and air mass flowrates given a certain power, equivalence ratio
# and NH3 mole fraction (a mixture of NH3 and H2 is assumed)
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

def CalcLHV(y_nh3):

    # This function calculates the LHV of the fuel given the NH3 mole fraction

    # Calculate the stoichiometric coefficients
    y_h2 = 1 - y_nh3
    y = y_h2*2 + y_nh3*3
    w = y_nh3
    x = 0.
    m = x + y/4                                 # O2/fuel molar basis

    # Calculate the fuel flowrate
    nH2O_s = y/2;                               # Stoichiometric moles of H2O produced
    DH_r   = -241.86*nH2O_s + 46.01*y_nh3;      # DH of global reaction kJ/mol
    LHV = -1000*DH_r/(y_h2*2 + y_nh3*17);       # LHV kJ/kg

    return LHV

# CALCULATE THE OUTLET TEMPERATURE OF THE NETWORK
def calc_outlet_temp(Power, y_nh3, phi_rich, phi_lean, T_in_rich, T_in_lean, T_amb=300.0, U=0.0, A=1.0):

    # This function calculates the outlet temperature of the network given:
    # Power: Power of the network (kW)
    # y_nh3: Mole fraction of NH3 in the fuel
    # phi_rich: Equivalence ratio of the rich combustion
    # phi_lean: Equivalence ratio of the lean combustion
    # T_in_rich: Inlet temperature of the rich combustion (K)
    # T_in_lean: Inlet temperature of the lean combustion (K)
    # T_amb: Ambient temperature (K)
    # U: Heat transfer coefficient (W/m2K)
    # A: Area of the reactor (m2)

    # Calculate the mass flowrate of the fuel and air
    m_fuel, m_air_lean = calc_inlet_mass(Power, y_nh3, phi_lean)
    m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)

    # Global mass of air
    m_air = m_air_rich + m_air_lean

    # Calculate the fraction of power reacted in the rich combustion
    y_h2 = 1 - y_nh3
    y = y_h2*2 + y_nh3*3
    w = y_nh3
    x = 0.
    m = x + y/4                                 # O2/fuel molar basis
    f = 0.79/0.21                               # Ratio between N2 and O2 in air
    af_mol  = m + m*f                           # Air/fuel stoichiometric molar basis
    af_mass = af_mol*29/(y_h2*2 + y_nh3*17)     # Air/fuel stoichiometric mass basis

    # In the rich zone, there is less air than the stoichiometric value. Therefore,
    # the fraction of fuel reacted is equal to the ratio between the mass of air
    # and the stoichiometric mass of air
    f_rich = m_air_rich/af_mass

    # The power reacted in the rich zone is equal to:
    Power_rich = f_rich*Power # kW

    # Create Cantera objects so that the calculation is more accurate
    fuel = ct.Solution('/Users/matteosavarese/Desktop/Dottorato/Github/Dakota_Interface/MacOS/Stagni_NH3/chem.cti')
    fuel_comp = 'H2:' + str(y_h2) + ', NH3:' + str(y_nh3)
    fuel.TPX = 340.0, 101325.0, fuel_comp

    # Create air
    air = ct.Solution('/Users/matteosavarese/Desktop/Dottorato/Github/Dakota_Interface/MacOS/Stagni_NH3/chem.cti')
    air.TPX = T_in_rich, 101325.0, 'O2:0.21, N2:0.79'

    # Use quantity objects so that you can mix the fuel and air
    qfuel = ct.Quantity(fuel, mass=m_fuel)
    qair_rich  = ct.Quantity(air, mass=m_air_rich)

    # Combustion products
    prod_rich = qfuel + qair_rich

    # Calculate the chemical equilibrium
    prod_rich.equilibrate('HP')
    T_out_rich = prod_rich.T

    # Now mix the remaining air with the products of the rich combustion
    qair_lean  = ct.Quantity(air, mass=m_air_lean)

    # Second stage
    prod_lean = prod_rich + qair_lean
    prod_lean.equilibrate('HP')

    # Calculate the outlet temperature of the lean combustion
    T_out_lean = prod_lean.T

    return T_out_rich, T_out_lean, prod_rich, prod_lean

# This function directly calculates the phi_lean as a function of the outlet temperature T
def CalcPhiLean(Tout, *parameters):

    T, Power, y_nh3, phi_rich, T_in_rich, T_in_lean = parameters

    # Calculate the mass flowrate of the fuel and air
    m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)

    # Use cantera to find equilibrium mixture of rich stage
    fuel = ct.Solution('/Users/matteosavarese/Desktop/Dottorato/Github/Dakota_Interface/MacOS/Stagni_NH3/chem.cti')
    fuel_comp = 'NH3:1'
    fuel.TPX = 340.0, 101325.0*5, fuel_comp

    # Create air
    air = ct.Solution('/Users/matteosavarese/Desktop/Dottorato/Github/Dakota_Interface/MacOS/Stagni_NH3/chem.cti')
    air.TPX = T_in_rich, 101325.0*5, 'O2:0.21, N2:0.79'

    # Use quantity objects so that you can mix the fuel and air
    qfuel = ct.Quantity(fuel, mass=m_fuel)
    qair_rich  = ct.Quantity(air, mass=m_air_rich)

    # Combustion products
    prod_rich = qfuel + qair_rich

    # Calculate the chemical equilibrium
    prod_rich.equilibrate('HP')

    # Do a manual thermal balance to find the mass of air in the lean zone
    # Calculate the mass of fuel from prod_rich
    y_nh3_rich = prod_rich.Y[fuel.species_index('NH3')]
    y_h2_rich = prod_rich.Y[fuel.species_index('H2')]
    m_fuel_lean = prod_rich.mass * (y_h2_rich + y_nh3_rich)
    
    # Calculate the power left
    LHV_NH3 =    CalcLHV(y_nh3) * 1000                       # J/kg
    LHV_H2 =   CalcLHV(1 - y_nh3) * 1000                   # J/kg
    LHV = (y_h2_rich * LHV_H2 + y_nh3_rich * LHV_NH3)/(y_h2_rich + y_nh3_rich)        # J/kg
    Power_lean = m_fuel_lean * LHV                          # W

    # Create a new quantity object for the lean zone
    m_air_lean = (Power_lean - prod_rich.mass * qair_rich.cp_mass * (Tout - prod_rich.T))/(qair_rich.cp_mass*(Tout - T_in_lean))

    # Now get the equivalence ratio of the lean zone
    # Calculate the fraction of power reacted in the rich combustion
    y_h2 = 1 - y_nh3
    y = y_h2*2 + y_nh3*3
    w = y_nh3
    x = 0.
    m = x + y/4                                 # O2/fuel molar basis
    f = 0.79/0.21                               # Ratio between N2 and O2 in air
    af_mol  = m + m*f                           # Air/fuel stoichiometric molar basis
    af_mass = af_mol*29/(y_h2*2 + y_nh3*17)     # Air/fuel stoichiometric mass basis

    excess_lean = m_air_lean/(af_mass*m_fuel) # Equivalence ratio of the lean zone
    excess_global = (m_air_lean + m_air_rich)/(af_mass*m_fuel)

    phi_lean = 1/(1+excess_lean)
    phi_global = 1/(1+excess_global)
     

    return phi_lean, phi_global, m_air_rich, m_air_lean

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
# Check if phi_rich is in the input name list
# Phi lean should be calculated given the outlet temperature of the flame
print(names)
if 'phi_rich' in names and 'T_cstr_1' in names:
    phi_rich = float(values[names.index('phi_rich')])
    T_cstr_1 = float(values[names.index('T_cstr_1')])
else:
    print('Error: phi_rich or T_cstr_1 not found in input file')
    sys.exit(1)

# Calculate phi lean
T_in_lean = T_cstr_1    # K
parameters = (Tout, Power, y_nh3, phi_rich, T_in_rich, T_in_lean)
phi_lean, phi_global, m_air_rich, m_air_lean = CalcPhiLean(Tout, *parameters)

# Get mass flowrates
M = calc_mass_flowrates(P, y_nh3, phi_rich, phi_global)

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
                    line[line.index(subitem)] = str(np.sum(M[0,:]))
                elif subitem =='M1':
                    line[line.index(subitem)] = str(np.sum(M[1,:]))
                elif subitem =='M2':
                    line[line.index(subitem)] = str(np.sum(M[:,2]))
                elif subitem =='M3':
                    line[line.index(subitem)] = str(np.sum(M[:,3]))
                elif subitem =='M4':
                    line[line.index(subitem)] = str(np.sum(M[:,4]))
                elif subitem =='M5':
                    line[line.index(subitem)] = str(np.sum(M[:,5]))
                elif subitem =='M6':
                    line[line.index(subitem)] = str(np.sum(M[:,6]))
                elif subitem =='M7':
                    line[line.index(subitem)] = str(np.sum(M[:,7]))
        
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
# Modify the mass flowrates
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
                if subitem =='M0':
                    line[line.index(subitem)] = str(np.sum(M[0,:]))
                elif subitem =='M1':
                    line[line.index(subitem)] = str(np.sum(M[1,:]))
                elif subitem =='M2':
                    line[line.index(subitem)] = str(np.sum(M[:,2]))
                elif subitem =='M3':
                    line[line.index(subitem)] = str(np.sum(M[:,3]))
                elif subitem =='M4':
                    line[line.index(subitem)] = str(np.sum(M[:,4]))
                elif subitem =='M5':
                    line[line.index(subitem)] = str(np.sum(M[:,5]))
                elif subitem =='M6':
                    line[line.index(subitem)] = str(np.sum(M[:,6]))
                elif subitem =='M7':
                    line[line.index(subitem)] = str(np.sum(M[:,7]))
            
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

    H2O_out = dt['H2O_x(18)'].values[0]
    O2_out = dt['O2_x(13)'].values[0]
    
    if s[0] == 'T':
        output_values.append(dt['T[K](5)'].values[-1])
    elif s[0] == 'NO':
        output_values.append(dt['NO_x(10)'].values[-1]*1e6)
    elif s[0] == 'NH3':
        output_values.append(dt['NH3_x(11)'].values[-1]*1e6)

if len(output_values) != len(output_names):
    print('Error: output_values and output_names have different lengths')
    sys.exit(1)

# Write outputs to file
output_dakota = sys.argv[2]
with open(output_dakota, 'w') as f:
    for i in range(len(output_values)):
        f.write(str(output_values[i]) + ' ' + output_names[i] + '\n')


