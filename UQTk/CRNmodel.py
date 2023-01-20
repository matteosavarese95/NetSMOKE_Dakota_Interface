#!/usr/bin/env python

# This is the black-box model that will be run by uq_pc.py as model_bb

# WARNING: 
# input parameters file is presumed to contains columns of T_CSTR_1, TAU_CSTR_2, H_CSTR_2, L_PFR
# x-cond file is presumed to contains columns of phi_rich, T_CSTR_3


import argparse
import os
import sys
from genericpath import isfile
import numpy as np
import math
import pandas as pd


# ---------------- USER DEFINED PARAMETERS ---------------- #
Pwr = 48.0        # kW
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

#def calc_phi_lean()
 


def CRN_fwd_model(upars, xcond, outfile=None):
    
    Nsamp = upars.shape[0]
    Npars = upars.shape[1]
    print("Number of samples found: %d" %Nsamp)
    print("Uncertain parameters dimension: %d" %Npars)
    
    Ncond = xcond.shape[0]
    Nx = xcond.shape[1]
    print("Number of x-conditions found: %d" %Ncond)
    print("x-conditions dimension: %d" %Nx)

    # Get the number of reactors and the list of files to open to replace strings
    print('Checking NetSmoke files to be open...\n')

    basefolder = 'CRN_model/'

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

	# WARNING: 
	# input parameters file is presumed to contains columns of T_CSTR_1, TAU_CSTR_2, H_CSTR_2, L_PFR
	# x-cond file is presumed to contains columns of phi_rich, T_CSTR_3
    # outputs will be 'T_cstr_2' ,'T_pfr_4', 'NO_pfr_4', 'NH3_pfr_4'
    
    parnames = ['phi_rich', 'T_CSTR_3','phi_lean','T_CSTR_1','TAU_CSTR_2','H_CSTR_2','L_PFR']
    
    response = []
    
    for s in range(Nsamp):
        #from upars
        T_CSTR_1 = upars[s,0]            
        TAU_CSTR_2 = upars[s,1]
        H_CSTR_2 = upars[s,2]
        L_PFR = upars[s,3]
    	
        
        all_outputs_one_samp = []

        for c in range(Ncond):
            #from x-cond:
            phi_rich = xcond[c,0]
            T_CSTR_3 = xcond[c,1]
            phi_lean =  0.4 #to be computed from phi_rich and phi_global
                        
            parvalues = [phi_rich,T_CSTR_3,phi_lean,T_CSTR_1,TAU_CSTR_2,H_CSTR_2,L_PFR]
            print("Running condition %d, sample %d" %(c,s))

            # Get mass flowrates
            M = calc_mass_flowrates(Pwr, y_nh3, phi_rich, phi_lean)
    		
            plh = []
            for i in range(5):
                for j in range(5):
                    plh.append('M'+str(i)+str(j))
    		
            #build folder structure for present case
            folder = 'CRN_c'+str(c)+'_s'+str(s)+'/'	

            os.system('mkdir '+folder )

            os.system('cp -r '+basefolder+'/* '+folder+'.')


    		# Modify each dictionary
            placeholders = ['M0','M1','M2','M3','M4'] + parnames
            newvals = [M[0,2], M[1,2], M[2,4], M[3,4], M[2,4]+M[3,4]] + parvalues
            for fname in fname_list:
                for p,v in zip(placeholders,newvals):
                    cmd = "sed -i '' -e 's/"+str(p)+"/"+str(v)+"/g' "+folder+fname
                    os.system(cmd)

            print('New inputs are ready')
            print('Now modifying the internal connections of the network...')
            
            # ---------------- WRITING THE MAIN INPUT.DIC ---------------- #
            
            # Modify main input dictionary
            
            fname = 'input.dic'
            for p,v in zip(plh+placeholders,M.reshape((25)).tolist()+newvals):
                cmd = "sed -i '' -e 's/"+str(p)+"/"+str(v)+"/g' "+folder+fname
                os.system(cmd)

            print('Main input file is ready. Run NetSMOKE')
            
            # ---------------- RUNNING THE SIMULATION ---------------- #
            os.chdir(mydir+'/'+folder)
            os.system('/Users/imacremote/Distributions/NetSmoke_Linux-Master/SeReNetSMOKEpp-master/bin/SeReNetSMOKEpp.sh --input input.dic')
            os.chdir(mydir)
            
            # ---------------- READING CRN OUTPUTS ---------------- #
            print('Reading outputs...')
            
            output_names = ['T_cstr_2' ,'T_pfr_4', 'NO_pfr_4', 'NH3_pfr_4']
            output_values = []
            for item in output_names:
                spl = item.split('_')
                outname = folder+'Output/Reactor.' + spl[2] + '/Output.out'
                dt = pd.read_csv(outname, sep = '\s+')
            
                H2O_out = dt['H2O_x(17)'].values[0]
                O2_out = dt['O2_x(15)'].values[0]
                
                if spl[0] == 'T':
                    output_values.append(dt['T[K](5)'].values[-1])
                elif spl[0] == 'NO':
                    output_values.append(dt['NO_x(21)'].values[-1]*1e6)
                elif spl[0] == 'NH3':
                    output_values.append(dt['NH3_x(32)'].values[-1]*1e6)
            
            all_outputs_one_samp.append(output_values)  #will have 
            
            # ---------------- WRITING OUTPUTS TO FILE ---------------- #
            if outfile:               
                with open(outfile, 'a') as f:
                    for i in range(len(output_values)):
                        f.write(str(output_values[i]) + '  ')
        
        # response is a Nsamp-rows and (Ncond*Noutputs)-columns matrix                
        response.append([item for sublist in all_outputs_one_samp for item in sublist]) #append flattened out list
        
        if outfile:
            with open(outfile, 'a') as f:
                f.write("\n")  
                
    return(response)
        
   

        
def main(argv):

    ## Parse input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outs', type=int,nargs='*',help="Range of indices of requested outputs (count from 1)")

    parser.add_argument("-i", "--input",   dest="input_parameters_file",   type=str,   default='ptrain.dat', help="Input parameters file, e.g. ptrain.dat")
    parser.add_argument("-x", "--xcond",   dest="x_conditions_file",       type=str,   default='xcond.dat',  help="X-conditions file, e.g. xcond.dat")
    parser.add_argument("-o", "--output",  dest="outputs_file",            type=str,   default='ytrain.dat', help="Outputs parameters file, e.g. ytrain.dat")
    args = parser.parse_args()

    if (os.path.isfile(args.input_parameters_file)):
        upars=np.loadtxt(args.input_parameters_file)  #(Nsamp x Npars)
    else:
        print('Error: %s not found' %args.input_parameters_file)
        sys.exit(1)
        

    if (os.path.isfile(args.x_conditions_file)):
        xcond=np.loadtxt(args.x_conditions_file)  #(Ncond x Nx)
    else:
        print('Error: %s not found' %args.x_conditions_file)
        sys.exit(1)

    
    if (os.path.isfile(args.outputs_file)):
        print('%s file already existing. New results will be appended to this file' %args.outputs_file)


    # CALL FORWARD MODEL
    
    response = CRN_fwd_model(upars, xcond, outfile=args.outputs_file)
            

if __name__ == "__main__":
   main(sys.argv[1:])


