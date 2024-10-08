#!/bin/bash -e
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
#=====================================================================================

# Script-example for 2d model ensemble generation

# Directory forward UQ scripts
UQPC=${UQTK_INS}/examples/uqpc
UQBIN=${UQTK_INS}/bin

# Create parameter ranges file
echo "600 950" > prange.dat
echo "0.04 0.06" >> prange.dat
echo "0.9411768 1.4117652" >> prange.dat
echo "2000 6000" >> prange.dat

NTRAIN=2048   # N sampling points
NVAL=128     # V sampling points
SAMPLING=rand  #or quad


#Generate x-conditions in xcond.dat

#TBD
#in python would be
#for i in np.arange(1,1.6,0.1):
#    for j in np.arange(800,1050,50):
#        print(i,j)



# Prepare random samples
${UQPC}/uq_pc.py -r offline_prep -p prange.dat -s $SAMPLING -n $NTRAIN -v $NVAL
#${UQPC}/uq_pc.py -r offline_prep -p prange.dat -s quad -n 4 

###################################################################################
###################################################################################


# Generate forward model evaluations 
python CRNmodel.py -i ptrain.dat -q qtrain.dat -o ytrain.dat -x xcond.dat 
python CRNmodel.py -i pval.dat   -q qval.dat   -o yval.dat   -x xcond.dat 
