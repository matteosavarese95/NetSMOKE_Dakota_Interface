# NetSMOKE_Dakota_Interface
This folder contains the MacOS (new version) and the Windows (old version) of the Dakota-NetSMOKE interface

## Python packages used for MacOS version
Python 3.10.8
pandas 1.1.3
numpy 1.23.1

## Dakota version for MacOS
dakota 6.16

## Instructions
### Dakota Input
dakota.in file: You should modify this file according to your preferences for Dakota settings.
In the section 'variables', under the voice 'descriptors', you specify the name of the parameters to be changed in the simulation.
No specific format is required, but the corresponding descriptors need to be present in the Netsmoke input files
e.g. if you specify:
'descriptor T_cstr_2'

in the Netsmoke input file (e.g. input.cstr.2.dic) you should replace:
@Temperature YYY.YY K; with
@Temperature T_cstr_2 K;

Note that the descriptors 'phi_rich' and 'phi_lean' should always be present in the dakota.in

Under the section 'responses', specify the desired descriptors for the output.
Please name those descriptors as:
T_cstr_2, NO_cstr_4 or similar, otherwise they will not be read by the interface



