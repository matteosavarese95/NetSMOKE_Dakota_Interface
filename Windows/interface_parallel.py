#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 10:17:49 2020

@author: andrea
"""
#import numpy as np
#import pandas as pd
#import sys
import os
import xmltodict as xtd
import copy
#from collections import ordereddict
import math
import re
import numpy as np

# ------------
# Preprocess to read network mass flowrates
#
# Mij is the n x n matrix with the net mass flowrate between each reactor
# the element ij is the mass flowrate from reactor i to reactor j

# We want to create a matrix alpha(ij) where the element ij is the fraction of
# mass from reactor i to reactor j
def get_mass_split(Mij, m_ext):
    
    # Check the dimensions
    (n,m) = np.shape(Mij)
    if n != m:
        print("Mij dimensions do not agree")
   
    # Initialize alpha
    alpha = np.copy(Mij)
    
    for i in range(n):
        # Check wether reactor i is an outlet reactor
        if m_ext[i,1] != 0:
            alpha[i,:] = Mij[i,:]/(sum(Mij[i,:]) + m_ext[i,1])
        
        else:
            alpha[i,:] = Mij[i,:]/sum(Mij[i,:])

        
    return alpha

# This function will check the mass balance of the network. m_ext is a n x 2
# array containing the external input in column 1 and the external output in column 2
def check_mass_balance(Mij, m_ext):
    
    # Check the dimensions
    (n,m) = np.shape(Mij)
    if n != m:
        print("Mij dimensions do not agree")
        
    mb = True
        
    # Check relative local mass imbalance
    li = np.zeros((n,1))
    for i in range(n):
        M_in  = sum(Mij[:,i])
        M_out = sum(Mij[i,:])
        M_in_tot = M_in + m_ext[i,0]
        M_out_tot = M_out + m_ext[i,1]
        li[i] = np.abs((M_in_tot - M_out_tot)/min(M_in_tot, M_out_tot))
        
        # Check the mass balance
        if li[i] > 1e-2:
            print("Local mass imbalance in reactor ", i, " > 1%")
            mb = False
            
    if mb == True:
        print('Mass balance within the tolerance')
    else:
        print('Check mass balance')
            
        
    return li, mb
        

# This function will get back the mass flowrates from the splitting matrix, 
# given the external inputs and outputs m_ext
from scipy import linalg
def get_mass_flowrates(alpha, m_ext):
    
    (n,m) = np.shape(alpha)
    if n != m:
        print("Mij dimensions do not agree")
        
    A = np.eye(n) - np.transpose(alpha)
    b = m_ext[:,0]
    
    # Mass flowrate from each reactors
    Msol = linalg.solve(A, b)
    
    Mij = np.copy(alpha)
    for i in range(n):
        Mij[i,:] = alpha[i,:]*Msol[i]
    
    # Check mass balance
    li, mb = check_mass_balance(Mij, m_ext)
    if mb == True:
        print('Mass flowrates succesfully found')
    else:
        print('Mass flowrates found do not satisfy mass balance')
        
    return Msol, Mij

def change_alpha(alpha, m_ext, ind_alpha_change):
    
    # Inputs:
        # alpha: matrix of the split ratio, modified by dakota
        # m_ext: array nr x 2 of the external inlets (1) and outlets (2)
        # ind_alpha_change = list of tuple of the indexes of the modified alpha
        
    alpha_new = np.copy(alpha)
    
    # Sum along the columns, should be 1 except for outlet reactors
    alpha_sum = np.sum(alpha, 1)
    
    # Now check if the sum is not one
    for i in range(len(alpha_sum)):
        
        # If the sum is not one and the reactor is not an outlet the untouched alpha
        # should be modified so that the sum is 1
        if alpha_sum[i] != 1 and m_ext[i,1] == 0.0:
            
            # Check the non-zero elements of alpha[i,:]
            non_zero_alpha = []
            for j in range(len(alpha_sum)):
                if alpha[i,j] != 0.0:
                    non_zero_alpha.append((i,j))
                    
            
            # Non zero alpha is a list of tuple with the non-zero elements of the
            # i-th row of alpha
            
            # Now check which elements of non_zero_alpha have been changed or not from
            # dakota by searching in ind_alpha_change
            ind_rows_c = []
            ind_rows_nc = []
            for j in range(len(non_zero_alpha)):
                if non_zero_alpha[j] in ind_alpha_change:
                    ind_rows_c.append(non_zero_alpha[j])
                
                elif non_zero_alpha[j] not in ind_alpha_change:
                    ind_rows_nc.append(non_zero_alpha[j])
                    
                    
            # Now check the dimension of ind_rows_nc
            if len(ind_rows_nc) == 1:
                
               # If the length of ind_rows_c is one:
                   if len(ind_rows_c) == 1:
                       alpha_new[ind_rows_nc[0]] = 1.0 - alpha[ind_rows_c[0]]
                       
                   else:
                       s = 0.0
                       for j in range(len(ind_rows_c)):
                           s = s + alpha(ind_rows_c[j])
                       
                       alpha_new[ind_rows_nc[0]] = 1.0 - s
                       
           
    
    return alpha_new

def calc_num_reactors(d):
    
    with open(r'{}'.format(d) + '\\toy_network_2\\input.dic') as dic:
        c = dic.readlines()
        index = [i for i, s in enumerate(c) if '@InternalConnections' in s]
        ind = index[0]
        r = []
        for s in range(1000):
            ind = ind + 1
            
            if len(c[ind].split('\t')) != 3:
                break
            else:
                string = c[ind].split('\t')
                r.append(int(string[0]))
                r.append(int(string[1]))
        
        nr = max(r) + 1
    return nr

                    
                
        

#--------------
#preprocessing
#--------------
def double_split(text_to_split):
    '''
    Useful function to completely clear array values extracted as string from dictionaries
    return:
        -values, strings array
    '''
    values = []
    split1 = text_to_split.split("\n")
    for s in split1:
        split2 = s.split(" ")
        for v in split2:
            values.append(v)
    return values

def read_input_files(d, nr):

    '''
    This function read files from input folder
    return:
        - kinetics_dict, dictionary parsed from xml file
        - r_names_dict,  dictionary parsed from xml file
        - factors_input, list of factors (strings)
        - reactions_input, list of reaction numbers (strings)
    '''

    with open(r'{}\..\kinetics'.format(d) + '\kinetics.xml') as kinetic_mechanism:
        km = kinetic_mechanism
        kinetics_dict = xtd.parse(km.read())


    with open(r'{}\..'.format(d) + '\\template.dic', mode="r") as os_t:
        template = os_t.readlines()


    # with open('input.initial.dic', mode="r") as os_t:
    #  template = os_t.readlines()
    #print(template)

    with open(r'{}'.format(d) + '\dakota_params.in') as f:
        # read first line
        a = f.readline().strip().split(" ")
        # extract number of parameters
        var_num = int(a[0])
        # initialize targets labels and suggested values
        p_values = [0] * var_num
        parameters = [0] * var_num
        # create a range for a for loop
        x = range(var_num)
        for n in x:
            # read line into a list
            b = f.readline().strip().split(" ")
            #assign value and label
            p_values[n] = float(b[0])
            parameters[n] = b[1]
            
    with open(r'{}'.format(d) + '\\toy_network_2\\input.dic') as dic:
        c = dic.readlines()
        index = [i for i, s in enumerate(c) if '@InternalConnections' in s]
        ind_ret = index[0]
        ind = index[0]
        M_net = np.zeros((nr,nr))
        # Get mass flowrates
        for s in range(1000):
            ind = ind + 1
            
            if len(c[ind].split('\t')) != 3:
                break
            else:
                string = c[ind].split('\t')
                r1 = int(string[0])
                r2 = int(string[1])
                try:
                    mf = float(string[2])
                except:
                    mf = float(string[2].split(';')[0])
                    
                M_net[r1, r2] = mf
        
        # Get inlet mass flowrates
        index = [i for i, s in enumerate(c) if '@InputStreams' in s]
        ind = index[0]
        m_ext = np.zeros((nr, 2))
        for s in range(1000):
            ind = ind + 1
            
            if len(c[ind].split('\t')) != 2:
                break
            else:
                string = c[ind].split('\t')
                r1 = int(string[0])
                try:
                    mf = float(string[1])
                    print(mf)
                except:
                    mf = float(string[1].split(';')[0])
                    
                m_ext[r1,0] = mf
        
        # Get outlet mass flowrates
        index = [i for i, s in enumerate(c) if '@OutputStreams' in s]
        ind = index[0]
        for s in range(1000):
            ind = ind + 1
            
            if len(c[ind].split('\t')) != 2:
                break
            else:
                string = c[ind].split('\t')
                r1 = int(string[0])
                try:
                    mf = float(string[1])
                except:
                    mf = float(string[1].split(';')[0])
                    
                m_ext[r1,1] = mf
                
    return kinetics_dict, parameters, p_values, template, M_net, m_ext, ind_ret




def write_direct(label , index, new_doc, p_values):
    '''
    This function insert a random parameter.
    '''
    #CREATING lnA INCREMENTED XML
    if label[0] == "lnA":

        lnA_string = new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['lnA']
        lnA_values = double_split(lnA_string)

        lnA_values[int(label[1].split("R")[1])] = p_values[index]
        new_lnA_string = list_to_string(lnA_values)
        new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['lnA'] = ""
        new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['lnA'] = new_lnA_string

    #CREATING Beta INCREMENTED XML
    elif label[0] == "Beta":

        beta_string = new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['Beta']
        beta_values = double_split(beta_string)

        beta_values[int(label[1].split("R")[1])] = p_values[index]
        new_beta_string = list_to_string(beta_values)
        new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['Beta'] = ""
        new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['Beta'] = new_beta_string


    #CREATING E over R INCREMENTED XML
    #
    elif label[0] == 'EoR':

        E_over_R_string = new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['E_over_R']
        E_over_R_values = double_split(E_over_R_string)

        E_over_R_values[int(label[1].split("R")[1])] = p_values[index]
        new_E_over_R_string = list_to_string(E_over_R_values)
        new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['E_over_R'] = ""
        new_doc['opensmoke']['Kinetics']['KineticParameters']['Direct']['E_over_R'] = new_E_over_R_string

    return new_doc

def list_to_string(array_values):
    unique_string = ""
    for v in array_values:
        unique_string = unique_string + str(v) + "\n"
    return unique_string

def write_classic_plog_lnA(cmd, idx, new_doc, p_values, how_many_PLOG):
    '''
    This function manipulate a dictionary object in order to create three new xml files
    updated the values of an optimized extended plog reaction.
    '''
    # extracts the index of the reaction and find it's position into the list of extended plog reactions
    extplog_string = new_doc['opensmoke']['Kinetics']['PressureLog']
    reactions_string = extplog_string.split("\n")
    reactions_single_string = ''.join(reactions_string[1:])
    reactions_list = reactions_single_string.split(" ")
    idx_extplogr = reactions_list.index(cmd[1])
    # aquires the node content: a string
    extplog_parameters_string = new_doc['opensmoke']['Kinetics']['KineticParameters']['PressureLog']
    # transforms the string in a list
    extplog_parameters_list = extplog_parameters_string.split("\n")
    # keeps only the odd positions of the list, because they correspond to the parameters value
    extplog_parameters_list_p = extplog_parameters_list[1::2]
    # acquires the string related to the wanted reaction and splits it in a list
    reaction_of_interest_string = extplog_parameters_list_p[idx_extplogr]
    reaction_of_interest_list   = reaction_of_interest_string.split(" ")
    # Sometimes there's a space at the end, so we remove just in case
    if reaction_of_interest_list[-1] == " ":
        reaction_of_interest_list =reaction_of_interest_list[0:-1]
    if reaction_of_interest_list[-1] == "":
        reaction_of_interest_list =reaction_of_interest_list[0:-1]

    Beta = float(p_values[idx])
    # Counts how many Arrhenius expressions there are in the PLOG
    N = (len(reaction_of_interest_list)-2)/4
    # Loop over the expression and replace
    for i in range(0,int(N)):
        new_A = math.log10(float(reaction_of_interest_list[i*4+1])) + Beta
        reaction_of_interest_list[i*4+1] = str(pow(10,new_A))


    # it recomposed the string from the list
    new_reaction_of_interest_string = " ".join(reaction_of_interest_list)
    # replace the string in the original list of extended plog reactions
    extplog_parameters_list[2*idx_extplogr+1] = new_reaction_of_interest_string
    # recompose original string
    new_extplog_parameters_string = "\n".join(extplog_parameters_list)
    #replace within the json file
    new_doc['opensmoke']['Kinetics']['KineticParameters']['PressureLog'] = new_extplog_parameters_string

    return new_doc

def update_input_conditions(new_doc, subs, k):
        for idx, row in enumerate(new_doc):
                a = re.sub(k, subs, row)
                new_doc[idx] = a
        return new_doc

#--------------
#preprocessing
#--------------
current_directory = os.path.dirname(os.path.abspath(__file__))

# Read the input files
nr = calc_num_reactors(current_directory)
kinetics_dict, parameters, p_values, template, M_net, m_ext, ind_net = read_input_files(current_directory, nr)
new_template = template
how_many_classic_PLOG = 0

print(parameters)

for idx, inp in enumerate(parameters):

        cmd = inp.split("_")

        if cmd[0] == "PLOG":

            how_many_classic_PLOG = how_many_classic_PLOG +1

print('Number of classic PLOGs = ' + str(how_many_classic_PLOG))

# Get split ratio
alpha = get_mass_split(M_net, m_ext)

counter = 0;

# Tuple list of index of alpha changed
ind_alpha_change = []

for idx, inp in enumerate(parameters):

    cmd = inp.split("_")

    if cmd[0]=="lnA":
        new_doc = write_direct(cmd , idx, kinetics_dict, p_values)
        kinetics_dict = new_doc
    elif cmd[0]=="EoR":
        new_doc = write_direct(cmd , idx, kinetics_dict, p_values)
        kinetics_dict = new_doc
    elif cmd[0]=="PLOG":
        new_doc = write_classic_plog_lnA(cmd , idx, kinetics_dict, p_values, how_many_classic_PLOG)
        kinetics_dict = new_doc
    elif cmd[0]=="SPLIT":
        new_key = "COMPLEMENT_" + cmd[2]
        document = copy.deepcopy(new_template)
        new_template = update_input_conditions(document, str(p_values[idx]), inp)
        document = copy.deepcopy(new_template)
        new_template = update_input_conditions(document, str(1-p_values[idx]), new_key)
        
    elif cmd[0] == 'V':
        rtype = cmd[1]
        rnum  = cmd[2]
        dic_name = 'input.' + rtype + '.' + rnum + '.dic'
             
        with open(current_directory + '\\toy_network_2\\' + dic_name, mode="r") as os_t:
            pfr_temp = os_t.readlines()
             
        document = copy.deepcopy(pfr_temp)
        new_template = update_input_conditions(document, str(p_values[idx]), inp)
         
        # Write new input file
        with open(current_directory + '\\toy_network_2\\' + dic_name, mode="w") as inp:
            inp.writelines("".join(new_template))
        
    elif cmd[0] == 'D':
        rtype = cmd[1]
        rnum  = cmd[2]
        dic_name = 'input.' + rtype + '.' + rnum + '.dic'
        
        with open(current_directory + '\\toy_network_2\\' + dic_name, mode="r") as os_t:
            pfr_temp = os_t.readlines()
             
        document = copy.deepcopy(pfr_temp)
        new_template = update_input_conditions(document, str(p_values[idx]), inp)
         
        # Write new input file
        with open(current_directory + '\\toy_network_2\\' + dic_name, mode="w") as inp:
            inp.writelines("".join(new_template))
            
    elif cmd[0] == 'alpha':
        r1 = int(cmd[1])
        r2 = int(cmd[2])
        alpha[r1,r2] = p_values[idx]
        ind_alpha_change.append((r1,r2))
      
    
    else:
        document = copy.deepcopy(new_template)
        new_template = update_input_conditions(document, str(p_values[idx]), inp)
    
        
    # Modify the alpha to make the mass balance ok:
    alpha_new = change_alpha(alpha, m_ext, ind_alpha_change)
            
    # Re-calculate new mass flowrates
    M_usless, M_new = get_mass_flowrates(alpha_new, m_ext)
    new_connections = []
    for i in range(nr):
        for j in range(nr):
            
            if M_new[i,j] != 0:
                new_connections.append(str(i) + ' ' + str(j) + ' ' + str(M_new[i,j]))
    
    new_connections[-1] = str(new_connections[-1]) + ';'                                                             
    with open(current_directory + '\\toy_network_2\\' + 'input.dic', mode="r") as inp:
        f = inp.readlines()
        for i, idx in enumerate(new_connections):
            f[ind_net + 1 + i] = idx + '\n'

    with open(current_directory + '\\toy_network_2\\' + 'input.dic', mode="w") as inp:
        inp.writelines("".join(f))
            
        
    
         
if not os.path.exists(r'{}\kinetics'.format(current_directory)):
    os.mkdir(r'{}\kinetics'.format(current_directory))
else:
    os.system('rmdir ' + r'{}\kinetics'.format(current_directory))
    os.mkdir(r'{}\kinetics'.format(current_directory))

with open(r'{}\kinetics'.format(current_directory) + "\kinetics.xml", 'w') as EoR_kinetics:
    EoR_kinetics.write(xtd.unparse(kinetics_dict, pretty=True))

os.system('copy ' + r'{}\..\kinetics\reaction_names.xml '.format(current_directory) + r'{}\kinetics'.format(current_directory))


#if os.path.isfile(r'{}'.format(current_directory) + '/input.dic'):
 #   os.system('rm -r input.dic')

