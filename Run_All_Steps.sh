################################################################################
#
# Author:       Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Date:         01.08.2025
# Type:         Bash shell + Python3
#
# Description:
# Script to run all the process steps for Aeolus2 analysis of Bump evolution.
#
# Important:
# First, it is necessary to have completed all steps from the installation of Aeolus2.0.
#
# To run this script:
# $ bash Run_All_Steps.sh
#
################################################################################

echo 
echo 
echo "Begin Run_All_Steps.sh"
echo
echo



echo '#######################'
echo '##  Aeolus2.0_model  ##'
echo '#######################'


cd Aeolus2.0_model/

conda activate dedalus2

echo
echo
echo '# Running Experiments'
echo
echo

time mpiexec -n 8 python3 aeolus2main.py --moist_convection=0 --external_forcing=0 --output_folder=../Data/Aeolus2.0_Output_Dry_Control
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=1 --external_forcing=0 --output_folder=../Data/Aeolus2.0_Output_MC_Control

time mpiexec -n 8 python3 aeolus2main.py --moist_convection=0 --external_forcing=barotropic --external_forcing_epsilon=0.1 --output_folder=../Data/Aeolus2.0_Output_Dry_Barotropic_Weak 
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=0 --external_forcing=barotropic --external_forcing_epsilon=0.2 --output_folder=../Data/Aeolus2.0_Output_Dry_Barotropic_Strong 
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=0 --external_forcing=baroclinic --external_forcing_epsilon=0.1 --output_folder=../Data/Aeolus2.0_Output_Dry_Baroclinic_Weak 
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=0 --external_forcing=baroclinic --external_forcing_epsilon=0.2 --output_folder=../Data/Aeolus2.0_Output_Dry_Baroclinic_Strong 

time mpiexec -n 8 python3 aeolus2main.py --moist_convection=1 --external_forcing=barotropic --external_forcing_epsilon=0.1 --output_folder=../Data/Aeolus2.0_Output_MC_Barotropic_Weak 
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=1 --external_forcing=barotropic --external_forcing_epsilon=0.2 --output_folder=../Data/Aeolus2.0_Output_MC_Barotropic_Strong 
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=1 --external_forcing=baroclinic --external_forcing_epsilon=0.1 --output_folder=../Data/Aeolus2.0_Output_MC_Baroclinic_Weak 
time mpiexec -n 8 python3 aeolus2main.py --moist_convection=1 --external_forcing=baroclinic --external_forcing_epsilon=0.2 --output_folder=../Data/Aeolus2.0_Output_MC_Baroclinic_Strong 



echo
echo
echo '# Running Postprocessing'
echo
echo

python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_Dry_Control           --output_file=Aeolus2.0_Output_Dry_Control.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_MC_Control            --output_file=Aeolus2.0_Output_MC_Control.nc

python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_Dry_Barotropic_Weak   --output_file=Aeolus2.0_Output_Dry_Barotropic_Weak.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_Dry_Barotropic_Strong --output_file=Aeolus2.0_Output_Dry_Barotropic_Strong.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_Dry_Baroclinic_Weak   --output_file=Aeolus2.0_Output_Dry_Baroclinic_Weak.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_Dry_Baroclinic_Strong --output_file=Aeolus2.0_Output_Dry_Baroclinic_Strong.nc

python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_MC_Barotropic_Weak    --output_file=Aeolus2.0_Output_MC_Barotropic_Weak.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_MC_Barotropic_Strong  --output_file=Aeolus2.0_Output_MC_Barotropic_Strong.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_MC_Baroclinic_Weak    --output_file=Aeolus2.0_Output_MC_Baroclinic_Weak.nc
python aeolus2_postprocessing.py --output_folder=../Data/Aeolus2.0_Output_MC_Baroclinic_Strong  --output_file=Aeolus2.0_Output_MC_Baroclinic_Strong.nc




echo
echo '########################'
echo '##  Analysis_Scripts  ##'
echo '########################'
echo


cd ../Analysis_Scripts/

echo
echo '# Running Aeolus2_Data_Diff_Remap.py'
echo

python Aeolus2_Data_Diff_Remap.py


echo
echo '# Running Aeolus2_Plot_Mapas.py'
echo

python Aeolus2_Plot_Mapas.py


echo
echo '# Running Aeolus2_Plot_Hovmoller.py'
echo

python Aeolus2_Plot_Hovmoller.py


echo
echo '# Running Aeolus2_Plot_Series.py'
echo

python Aeolus2_Plot_Series.py










echo 
echo 
echo "End Run_All_Steps.sh"
echo
echo


