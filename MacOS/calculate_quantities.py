import sys
import os
from genericpath import isfile
import numpy as np
import pandas as pd
import cantera as ct

# ---------------- NEW MODEL STRUCTURE ---------------- #
# REACTOR 0: Fuel inlet (fake reactor), T = T_fuel (300 K) fixed
# REACTOR 1: Air inlet (fake reactor), T = T_rich  (800 - 1000 K) 
# REACTOR 2: Ignition flame PSR, adiabatic reactor, small volume (~30cm3)
# REACTOR 3: Quenching reactor, non-adiabatic, U = 10-20 W/m2/K, small volume (in theory should not be reactive)
# REACTOR 4: First post-flame PFR (only rich mixture) it simulates the duct (tau ~ 0.02s)
# REACTOR 5: First lean combustion reactor (tau ~ 0.001 s), rich mixture reacts with a portion of the lean mixture
# REACTOR 6: Second lean combustion reactor (stagnation zone), rich mixture reacts with a portion of the lean mixture
# REACTOR 7: Third lean combustion reactor, now is isothermal (T = 1200 K) because of the heat loss in the quenching

# ---------------- USER DEFINED PARAMETERS ---------------- #
P = 96.0                # kW
y_nh3 = 1.0             # mole fraction
Tout  = 1170.0          # K
Power = 96.0            # kW
T_in_rich = 300.0       # K (fuel temperature)

# ---------------- USER DEFINED PARAMETERS ---------------- #
# If UQ of phi_rich and T_cstr_1 is not carried out, then the following values are used
phi_rich = 1.0
T_cstr_1 = 1000.0

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

    Power, y_nh3, phi_rich, T_in_rich, T_in_lean = parameters

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

# ---------------- CALCULATE QUANTITIES ---------------- #
Power = 96.0
y_nh3 = 1.0
phi_rich = 1.0
phi_global = 0.172
T_in_rich = 800
T_in_lean = 800

# Calculate mass flowrates by assigning phi
m_fuel, m_air_rich = calc_inlet_mass(Power, y_nh3, phi_rich)
m_fuel, m_air_tot  = calc_inlet_mass(Power, y_nh3, phi_global)

print('m_fuel = ' + str(m_fuel) + ' kg/s')
print('m_air_rich = ' + str(m_air_rich) + ' kg/s')
print('m_air_tot = ' + str(m_air_tot) + ' kg/s')

parameters = (Power, y_nh3, phi_rich, T_in_rich, T_in_lean)
phi_lean, phi_global, m_air_rich, m_air_lean = CalcPhiLean(Tout, *parameters)
print('phi global = ' + str(phi_global))


