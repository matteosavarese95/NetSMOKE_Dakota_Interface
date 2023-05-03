#!/usr/bin/env python

# This is the black-box model that will be run by uq_pc.py as model_bb
# Alternatively, the function response=CRN_fwd_model(upars, xcond) is callable inline

# WARNING: 
# input parameters file is presumed to contains columns of U_CSTR_3, Fa1, Fa3
# x-cond file is presumed to contains columns of phi_rich, T_CSTR_1

# ---------------- NEW MODEL STRUCTURE ---------------- #
# REACTOR 0: Fuel inlet (fake reactor), T = T_fuel (300 K) fixed
# REACTOR 1: Air inlet (fake reactor), T = T_rich  (800 - 1000 K) 
# REACTOR 2: Ignition flame PSR, adiabatic reactor, small volume (~30cm3)
# REACTOR 3: Quenching reactor, non-adiabatic, U = 10-20 W/m2/K, small volume (in theory should not be reactive)
# REACTOR 4: First post-flame PFR (only rich mixture) it simulates the duct (tau ~ 0.02s)
# REACTOR 5: First lean combustion reactor (tau ~ 0.001 s), rich mixture reacts with a portion of the lean mixture
# REACTOR 6: Second lean combustion reactor (stagnation zone), rich mixture reacts with a portion of the lean mixture
# REACTOR 7: Third lean combustion reactor, now is isothermal (T = 1200 K) because of the heat loss in the quenching


import argparse
import os
import sys
from genericpath import isfile
import numpy as np
import math
import pandas as pd
import cantera as ct


basefolder = 'CRN_Model_2/'
path_to_netsmoke = '/Users/imacremote/Distributions/NetSmoke_Linux-Master/SeReNetSMOKEpp-master/bin'
#chemfile = 'Otomo/otomo.yaml'
chemfile = 'Stagni_NH3/chem.cti'

# ---------------- USER DEFINED PARAMETERS ---------------- #
Pwr = 96.0              # kW
y_nh3 = 1.0             # mole fraction
Tout = 1225.0           # Desired T out
T_in_rich = 300.0       # K (fuel temperature)
phi_rich = 1.0          # equivalence ratio



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
def calc_mass_flowrates(Power, y_nh3, phi_rich, phi_global, stagn1=0.1, stagn2=0.15, entr=0.976):

    # This function calculates the mass flow rates across the reactors 
    # and stores them in a Matrix M[i,j] where M[i,j] is the mass flow rate
    # from reactor i to reactor j. Please note that this function defines
    # the internal structure of the reactor network. The user should modify
    # this function to change the reactor network model.
    # stagn1, stagn2 are the fractions of the bypass air that goes to the
    # 1st stagnation reactor and 2nd stagnation reactor respectively.
    # They can be specified by the user from the Dakota input file.

    # Global quantities
    m_fuel, m_air_tot = calc_inlet_mass(Power, y_nh3, phi_global)
    m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)
    m_air_bypass = m_air_tot - m_air_rich

    # Initialize mass flowrates matrix
    M = np.zeros((8,8))

    # First combustion zone
    M[0,2] = m_fuel                 # To flame PSR
    M[1,2] = m_air_rich*entr        # To flame PSR (we consider the fraction of air that is entrained)

    # The quenching reactor is number 7
    M[2,3] = m_fuel + m_air_rich*entr       # From ignition to quenching PSR
    M[3,4] = m_fuel + m_air_rich*entr       # From quenching to FIRST POST-FLAME flame PFR
    M[4,5] = m_fuel + m_air_rich*entr       # From flame to first stagnation reactor

    # Stagnation zone
    stagn3 = 1 - stagn1 - stagn2
    M[1,5] = m_air_bypass*stagn1 + m_air_rich*(1.0-entr) # From inlet air to first stagnation reactor
    M[1,6] = m_air_bypass*stagn2                         # From inlet air to second stagnation reactor
    M[1,7] = m_air_bypass*stagn3                         # From inlet air to outlet reactor

    # Sequential mass flowrates
    M[5,6] = m_fuel + m_air_rich + m_air_bypass*stagn1
    M[6,7] = m_fuel + m_air_rich + m_air_bypass*stagn1 + m_air_bypass*stagn2

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
    fuel = ct.Solution(chemfile)
    fuel_comp = 'H2:' + str(y_h2) + ', NH3:' + str(y_nh3)
    fuel.TPX = 340.0, 101325.0, fuel_comp

    # Create air
    air = ct.Solution(chefile)
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
    fuel = ct.Solution(chemfile)
    fuel_comp = 'NH3:1'
    fuel.TPX = 300.0, 101325.0*5, fuel_comp

    # Create air
    air = ct.Solution(chemfile)
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
    LHV_H2  =    CalcLHV(1 - y_nh3) * 1000                   # J/kg
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

def CalcAirMassNew(Tout, Tair):

    # Properties of the fuel
    y_nh3 = 1.0
    LHV = CalcLHV(y_nh3)

    # Specific heat of the burned mixture
    cp = 1.173 # kJ/kg K

    # Calculate the mass of fuel and the mass of rich air
    Power = 96.0
    phi_rich = 1.0
    m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)

    # Define the residual function
    def res(m_air):

        Tmix = (m_fuel * 300 + m_air * Tair)/(m_fuel + m_air)
        mcalc = Power/(cp*(Tout - Tmix))

        return mcalc - m_air

    # Calculate the total mass 
    m_start = Power/(cp*(Tout - Tair))

    # Solve the equation
    from scipy.optimize import fsolve
    m_air_tot = fsolve(res, m_start)[0]

    return m_air_tot, m_air_rich

# Calculate the mass flowrates of the network
def calc_mass_flowrates_new(Power, y_nh3, phi_rich, Tout, Tair, stagn1=0.1, stagn2=0.15, entr=0.976):

    # This function calculates the mass flow rates across the reactors 
    # and stores them in a Matrix M[i,j] where M[i,j] is the mass flow rate
    # from reactor i to reactor j. Please note that this function defines
    # the internal structure of the reactor network. The user should modify
    # this function to change the reactor network model.
    # stagn1, stagn2 are the fractions of the bypass air that goes to the
    # 1st stagnation reactor and 2nd stagnation reactor respectively.
    # They can be specified by the user from the Dakota input file.

    # Calculate total mass of air
    m_air_tot, m_air_rich = CalcAirMassNew(Tout, Tair)

    # Global quantities
    m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)

    m_air_bypass = m_air_tot - m_air_rich

    # Initialize mass flowrates matrix
    M = np.zeros((8,8))

    # First combustion zone
    M[0,2] = m_fuel                 # To flame PSR
    M[1,2] = m_air_rich*entr        # To flame PSR (we consider the fraction of air that is entrained)

    # The quenching reactor is number 7
    M[2,3] = m_fuel + m_air_rich*entr       # From ignition to quenching PSR
    M[3,4] = m_fuel + m_air_rich*entr       # From quenching to FIRST POST-FLAME flame PFR
    M[4,5] = m_fuel + m_air_rich*entr       # From flame to first stagnation reactor

    # Stagnation zone
    stagn3 = 1 - stagn1 - stagn2
    M[1,5] = m_air_bypass*stagn1 + m_air_rich*(1.0-entr) # From inlet air to first stagnation reactor
    M[1,6] = m_air_bypass*stagn2                         # From inlet air to second stagnation reactor
    M[1,7] = m_air_bypass*stagn3                         # From inlet air to outlet reactor

    # Sequential mass flowrates
    M[5,6] = m_fuel + m_air_rich + m_air_bypass*stagn1
    M[6,7] = m_fuel + m_air_rich + m_air_bypass*stagn1 + m_air_bypass*stagn2

    return M




def CRN_fwd_model(upars, xcond, outfile=None):
    
    Nsamp = upars.shape[0]
    Npars = upars.shape[1]
    print("Number of samples found: %d" %Nsamp)
    print("Uncertain parameters dimension: %d" %Npars)
    
    Ncond = xcond.shape[0]
    Nx = xcond.shape[1]
    print("Number of x-conditions found: %d" %Ncond)
    print("x-conditions dimension: %d" %Nx)


    if (os.path.isfile('discarded_backdoor.dat')):
        os.system('mv discarded_backdoor.dat discarded_$(date +%Y-%m-%dT%H%M%S).dat')
    with open('discarded_backdoor.dat', 'w') as fb:
        fb.write("# discarded samples. Count starts from 0. \n")
    
    # Get the number of reactors and the list of files to open to replace strings
    print('Checking NetSmoke files to be open...\n')


    nr = 0
    fname_list = []
    for i in range(1000):
        if isfile(basefolder+'input.cstr.'+str(i)+'.dic'):
            fname_list.append('input.cstr.' + str(i) + '.dic')
            nr += 1
        elif isfile(basefolder+'input.pfr.'+str(i)+'.dic'):
            fname_list.append('input.pfr.' + str(i) + '.dic')
            nr += 1
        else:
            print('Found', nr, 'reactors')
            break
    
    # Check if main input is present
    if not isfile(basefolder+'input.dic'):
        print('Error: input.dic not found')
        sys.exit(1)
            
    mydir = os.getcwd()

    parnames = ['phi_rich','T_cstr_1','T_cstr_3','V_cstr_2','U_cstr_3', 'V_cstr_5', 'V_cstr_6', 'Fa1', 'Fa3', 'Fa4']
    
    response = []
    discarded = []
    
    for s in range(Nsamp):
        #from upars
        U_CSTR_3 = upars[s,0]
        Fa1 = upars[s,1]
        Fa3 = upars[s,2]
        
        #Fa3 = 0.1    #stagn1
        Fa4 = 0.15   #stagn2
        V_CSTR_2 = 75 
        V_CSTR_5 = 800
        V_CSTR_6 = 400
        #Fa1 = 0.976  #entr
    	
        
        all_outputs_one_samp = []

        for c in range(Ncond):
            #from x-cond:
            phi_rich = xcond[c,0]
            T_CSTR_1 = xcond[c,1]
            T_CSTR_3 = T_CSTR_1    
            
            T_in_lean = T_CSTR_1
            T_in_rich = T_in_lean
            #parameters = (Tout, Pwr, y_nh3, phi_rich, T_in_rich, T_in_lean)
            #phi_lean, phi_global, m_air_rich, m_air_lean = CalcPhiLean(Tout, *parameters)
            # Get mass flowrates
            M = calc_mass_flowrates_new(Pwr, y_nh3, phi_rich, Tout, T_CSTR_1, Fa3, Fa4, Fa1)


                    
            parvalues = [phi_rich,T_CSTR_1,T_CSTR_3,V_CSTR_2,U_CSTR_3,V_CSTR_5,V_CSTR_6,Fa1,Fa3,Fa4]
            print('####################################################################################')
            print("Running condition %d, sample %d" %(c,s))
            print(f'phi_rich {phi_rich}, T_in {T_CSTR_1}')
            print(f'U_CSTR_3 {U_CSTR_3}, Fa1 {Fa1}, Fa3 {Fa3} ')

    		
            plh = []
            for i in range(len(M)):
                for j in range(len(M)):
                    plh.append('M'+str(i)+str(j))
    		
            #build folder structure for present case
            folder = 'CRN_c'+str(c)+'_s'+str(s)+'/'	

            os.system('mkdir '+folder )

            os.system('cp -r '+basefolder+'/* '+folder+'.')


    		# Modify each dictionary
            placeholders = ['M0','M1','M2','M3','M4','M5','M6','M7'] + parnames
            newvals = [np.sum(M[0,:]), np.sum(M[1,:]), np.sum(M[:,2]), np.sum(M[:,3]), np.sum(M[:,4]), np.sum(M[:,5]), np.sum(M[:,6]), np.sum(M[:,7])] + parvalues
            for fname in fname_list:
                for p,v in zip(placeholders,newvals):
                    cmd = "sed -i '' -e 's/"+str(p)+"/"+str(v)+"/g' "+folder+fname
                    os.system(cmd)

            print('New inputs are ready')
            print('Now modifying the internal connections of the network...')
            
            # ---------------- WRITING THE MAIN INPUT.DIC ---------------- #
            
            # Modify main input dictionary
            
            fname = 'input.dic'
            for p,v in zip(plh+placeholders,M.reshape((M.size)).tolist()+newvals):
                cmd = "sed -i '' -e 's/"+str(p)+"/"+str(v)+"/g' "+folder+fname
                os.system(cmd)

            print('Main input file is ready. Run NetSMOKE')
            
            # ---------------- RUNNING THE SIMULATION ---------------- #
            os.chdir(mydir+'/'+folder)
            os.system(f'{path_to_netsmoke}/SeReNetSMOKEpp.sh --input input.dic > netsmoke.out')
            os.chdir(mydir)
            
            # ---------------- READING CRN OUTPUTS ---------------- #
            print('Reading outputs...')
            
            #output_names = ['T_cstr_2', 'T_cstr_6', 'NO_cstr_6', 'NH3_cstr_6']
            output_names = ['NO_pfr_7', 'NH3_pfr_7']
            output_values = []
            for item in output_names:
                spl = item.split('_')
                outname = folder+'Output/Reactor.' + spl[2] + '/Output.out'
                dt = pd.read_csv(outname, sep = '\s+', on_bad_lines='skip')
            
                
                if spl[0] == 'T':
                    output_values.append(dt['T[K](5)'].values[-1])
                elif spl[0] == 'NO':
                    output_values.append(dt['NO_x(10)'].values[-1]*1e6)  #ppm
                elif spl[0] == 'NH3':
                    output_values.append(dt['NH3_x(11)'].values[-1]*1e6)  #ppm
            
            all_outputs_one_samp.append(output_values)  
            
           
            
        # response is a Nsamp-rows and (Ncond*Noutputs)-columns matrix  
        
        #add this sample response, unless the output is not good
        #if (any(np.array(all_outputs_one_samp)[:,1] < Tout*0.95)):
        #    print(f'!!!!!!!!! Discarding sample {s} because Tout is too low')
        #    discarded.append(s)
        #    with open('discarded_backdoor.dat', 'a') as fb:
        #        fb.write(str(s))
        #        fb.write("\n")
        #else:
        all_outputs_one_samp_flat = [item for sublist in all_outputs_one_samp for item in sublist]            
        response.append(all_outputs_one_samp_flat) #append flattened out list
        
        # ---------------- WRITING OUTPUTS TO FILE ---------------- #
        if outfile:
            with open(outfile, 'a') as f:
                f.write(' '.join(str(x) for x in all_outputs_one_samp_flat))
                f.write("\n")  
                
    return(response, discarded)


        
def main(argv):

    ## Parse input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outs', type=int,nargs='*',help="Range of indices of requested outputs (count from 1)")

    parser.add_argument("-i", "--input",   dest="input_parameters_file",   type=str,   default='ptrain.dat', help="Input parameters file, e.g. ptrain.dat")
    parser.add_argument("-q", "--qsamp",   dest="input_germs_samples_file",type=str,   default='qtrain.dat', help="Input germ samples file, e.g. qtrain.dat")
    parser.add_argument("-x", "--xcond",   dest="x_conditions_file",       type=str,   default='xcond.dat',  help="X-conditions file, e.g. xcond.dat")
    parser.add_argument("-o", "--output",  dest="outputs_file",            type=str,   default='ytrain.dat', help="Outputs parameters file, e.g. ytrain.dat")
    args = parser.parse_args()

    if (os.path.isfile(args.input_parameters_file)):
        upars=np.loadtxt(args.input_parameters_file)  #(Nsamp x Npars)
    else:
        print('Error: %s not found' %args.input_parameters_file)
        sys.exit(1)
 
 
    if (os.path.isfile(args.input_germs_samples_file)):
        qpars=np.loadtxt(args.input_germs_samples_file)  #(Nsamp x Npars)
    else:
        print('Error: %s not found' %args.input_germs_samples_file)
        sys.exit(1)       


    if (os.path.isfile(args.x_conditions_file)):
        xcond=np.loadtxt(args.x_conditions_file)  #(Ncond x Nx)
    else:
        print('Error: %s not found' %args.x_conditions_file)
        sys.exit(1)


    # take care of previously generated samples
    if (os.path.isdir('CRN_c0_s0')):
        os.system('mkdir old_samples')
        os.system('mv CRN_c* old_samples')
        if (os.path.isfile(args.outputs_file)):
            os.system(f'mv {args.outputs_file} old_samples/.')
        os.system('mv old_samples old_samples_$(date +%Y-%m-%dT%H%M%S)')
        
        
    if (os.path.isfile(args.outputs_file)):
        os.system(f'mv {args.outputs_file} {args.outputs_file}_$(date +%Y-%m-%dT%H%M%S)')
        

    # CALL FORWARD MODEL
    
    response, discarded = CRN_fwd_model(upars, xcond, outfile=args.outputs_file)
    
    # take care of failed samples
    if len(discarded):
        print('updating input files removing failed samples ...')
        upars_survived = np.delete(upars, discarded, axis=0)
        qpars_survived = np.delete(qpars, discarded, axis=0)
        os.system(f'mv {args.input_parameters_file} {args.input_parameters_file}.original')
        os.system(f'mv {args.input_germs_samples_file} {args.input_germs_samples_file}.original')
        np.savetxt(args.input_parameters_file, upars_survived)
        np.savetxt(args.input_germs_samples_file, qpars_survived)

            

if __name__ == "__main__":
   main(sys.argv[1:])


