#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:35:48 2020

@author: andreac
""" 
import pandas as pd
import numpy as np
import sys
import os

#def interpolate_NO(vector_1, vector_2, vector_3, value):

 #   for idx, i in enumerate(vector_1):
  #      if i > value:
   #         pos_max = idx
    #        break

    #for idx, i in enumerate(vector_1):
     #   if i < value:
      #      pos_min = idx
            #break

    #result2 = vector_2[pos_min] + (value-vector_1[pos_min])*(vector_2[pos_max]-vector_2[pos_min])/(vector_1[pos_max]-vector_1[pos_min])
    #linearly interpolate the vector_2 to find response to value
    #result3 = vector_3[pos_min] + (value-vector_1[pos_min])*(vector_3[pos_max]-vector_3[pos_min])/(vector_1[pos_max]-vector_1[pos_min])
    #linearly interpolate the vector_2 to find response to value
    #result = result2/(1-result3)

    #return result
    

current_directory = os.path.dirname(os.path.abspath(__file__))
output_file_name  = sys.argv[1]

#dt = pd.read_csv('os_out/Solution.final.out', sep = '\s+', header = 'infer', dtype = None, skipinitialspace=True)
#dt = pd.read_csv('os_out/Solution.final.out', sep = '\s+', usecols = [1, 285, 26], names = ["x[cm](2)", "NO_x(286)", "H2O_x(27)"], header = 0, dtype = None, skipinitialspace = False)
#dt = pd.read_csv('toy_network/Output/Reactor.1/Output.out', sep = '\s+', usecols = [48], names = ["NO_x(49)"], header = 0, dtype = None, skipinitialspace = False)
dt = pd.read_csv('toy_network_2/Output/Reactor.3/Output.out', sep = '\s+')

labels_results = ['T[K](5)', 'NO_x(21)', 'NH3_x(32)']
results = [None]*len(labels_results)

H2O_out = dt['H2O_x(17)'].values[0]
O2_out = dt['O2_x(15)'].values[0]

for idx,i in enumerate(labels_results):
    results[idx] = dt[i].values[0]
    
    if labels_results[idx] == 'NO_x(21)':
        results[idx] = 1e6 * results[idx]/(1-H2O_out) * (21 - 3)/(21 - O2_out*100)
    elif labels_results[idx] == 'NH3_x(32)':
        results[idx] = 1e6 * results[idx]/(1-H2O_out) * (21 - 3)/(21 - O2_out*100)
        
    #print(np.argmin(abs(dist_burner[pos]-i)))
print(results)

output_names = ['T_out', 'NO_out', 'NH3_out']

fh = open(r'{}/../'.format(current_directory) + output_file_name, "w")
for idx,i in enumerate(results):
    fh.write(repr(i) + '\t' + output_names[idx] + '\n')  
fh.close()
