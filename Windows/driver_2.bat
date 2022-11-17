::!/bin/bash

:: $1 and $2 are special variables in bash that contain the 1st and 2nd
:: command line arguments to the script, which are the names of the
:: Dakota parameters and results files, respectively.
:: fork interface runs: ./driver.sh params.in.n results.out.n
:: num=$(echo %1 | awk -F. '{print $NF}')         

:: Extract only the last 7 characters
SET string=%1
::SET num=%string:~-1%
set "num=%string:.=" & set "num=%"

:: assign current directoty to top tree directory
SET topdir=%cd%

:: create a workdir name for the specific thread
::workdir=%topdir\workdir.%num
SET workdir=%topdir%\workdir.%num%

:: make dir
md %workdir%

:: copy the params.in.n to the workdir
copy %topdir%\%string% %workdir%\dakota_params.in
copy interface_parallel.py %workdir%\interface_parallel.py
copy pp_pro_flame.py %workdir%\pp_pro_flame.py
md %workdir%\toy_network_2
xcopy toy_network_2 %workdir%\toy_network_2 

:: change directory to thread directory
cd %workdir%

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::
:::: Pre-processing Phase -- Generate/configure an input file for your simulation
::::  by substiting in parameter values from the Dakota paramters file.
::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

python interface_parallel.py
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::
:::: Execution Phase -- Run your simulation
::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: Enter CRN folder
cd toy_network_2

:: Run netsmoke
"C:\Users\matte\Desktop\Dottorato\NetSmoke\exe\SeReNetSMOKE++.exe" --input input.dic

:: Come back to workdir
cd ..\
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::
:::: Post-processing Phase -- Extract (or calculate) quantities of interest
::::  from your simulation's output and write them to a properly-formatted
::::  Dakota results file.
::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
python pp_pro_flame.py %2

cd ..\
