''' 
This is the configuration file for Aeolus2. 
Here, key parameters and user settings are defined. 
This file is then read by the main script, aeolus2main.py.

Structure:

    Path specifications: defines where to read/write data and how to name output files.
    User settings: submodule activations.
    Scaling and physical parameters.
    Discretization parameters.
    Integration parameters.
'''


### Imports
import numpy  as np
import math   as ma
import glob
import os


### path specification ###

output_folder            = '../Data/Aeolus2.0_Output_Dry_Control'				# (moist_convection = 0, external_forcing = 0)
# output_folder          = '../Data/Aeolus2.0_Output_Dry_Barotropic_Weak'		# (moist_convection = 0, external_forcing = 'barotropic', external_forcing_epsilon = 0.1)
# output_folder          = '../Data/Aeolus2.0_Output_Dry_Barotropic_Strong'		# (moist_convection = 0, external_forcing = 'barotropic', external_forcing_epsilon = 0.2)
# output_folder          = '../Data/Aeolus2.0_Output_Dry_Baroclinic_Weak'		# (moist_convection = 0, external_forcing = 'baroclinic', external_forcing_epsilon = 0.1)
# output_folder          = '../Data/Aeolus2.0_Output_Dry_Baroclinic_Strong'		# (moist_convection = 0, external_forcing = 'baroclinic', external_forcing_epsilon = 0.2)
# output_folder          = '../Data/Aeolus2.0_Output_MC_Control'				# (moist_convection = 1, external_forcing = 0)
# output_folder          = '../Data/Aeolus2.0_Output_MC_Barotropic_Weak'		# (moist_convection = 1, external_forcing = 'barotropic', external_forcing_epsilon = 0.1)
# output_folder          = '../Data/Aeolus2.0_Output_MC_Barotropic_Strong'		# (moist_convection = 1, external_forcing = 'barotropic', external_forcing_epsilon = 0.2)
# output_folder          = '../Data/Aeolus2.0_Output_MC_Baroclinic_Weak'		# (moist_convection = 1, external_forcing = 'baroclinic', external_forcing_epsilon = 0.1)
# output_folder          = '../Data/Aeolus2.0_Output_MC_Baroclinic_Strong'		# (moist_convection = 1, external_forcing = 'baroclinic', external_forcing_epsilon = 0.2)

output_file              = 'output_'     				# output_ --> [output_00001.npz, ..., output_00124.npz]

restart_from_output      = 0							# 1: To restart from the last file (Ex. output/output_0000N.npz) | 0: to use Obs path 

n_iterations             = 23200            	  		# total number of iterations (Ex: 7000 -> 30.5432 [days] = 35.4782 [nondim. days] // 23200 = Time 100.976 [days] = 117.2912 [nondim. days] -> 408 files)
cadence                  = 4.0                   		# number of output files per day 


#    --------------------- *** User setting & parameterization *** -----------
smooth_run               = 1                     		# 1: latitudinal fine resolution (768*384) with 1 min numerical time step; 
fast_run                 = 0                     		# 1: time step=step min, some second oreder nonliear terms have been modified, resolution: (384*192)
super_fast               = 0                     		# 1: time step=step min, some second oreder nonliear terms have been modified, resolution: (192*96)   

moist_convection         = 0
topography               = 1                     		# 1: inclusion of topography
restart_run              = 1                     		# 1: restart from previous run except TQG_adjustment
Radiative_Transfer       = 1                     		# 1: RRTMG_LW/SW: rapid radiative transfer model that utilizes the correlated-k approach to calculate longwave fluxes
unparallel_TQG_adj       = 0                     		# 1: provides fully nonlinear thermo-quasi-geostrophic balanced state of initial state or any other initialization
													    # In this case radiative transfer, and moist-convection are off.
													
external_forcing         = 0		            		# Any arbitrary forcing in the main script // Updated options: 0 (no forcing), 'baroclinic' (lower layer forcing), 'barotropic' (lower and upper layers forcing)
external_forcing_epsilon = 0                     		# maximum amplitude for the artificial buoyance anomaly forcing. E.g.: 0.1, 0.2, 0.3 
                      
Newtonian_cooling        = 1                     		# 1: Newtonian cooling is active, 0: thermal relaxation of bh toward initial BH level
summer_solstice          = 0                     		# 1: just insolation of summer solstice will be imposed
winter_solstice          = 0                     		# 1: just insolation of winter solstice will be imposed  

# Set the type of insolation forcing by choosing one of the following flags or defining your own type of insolation:
Insolation_Exclusive     = 1
Anomaly_Mean_Critical    = 0                     		# 1 indicates an insolation anomaly significantly deviating from the mean value.
#Anomaly_Ref_Critical    = 0                    		# 1 indicates an insolation anomaly significantly deviating from the reference level.
# ------------------------------------------------------------------

spin_up_days             = 5.0                   		# Sets the length of the spin-up period in days required for the model to initialize and stabilize before the main simulation. 
CO2ppmv                  = 280. 
T_s0                     = 250.0277              		# 262.1009#280.
P_s                      = 1000.0                		# [hPa] 
H0_dim                   = 10000.0               		# Needs to be adapted by corresponding scaling. Scale of (pseudo)-height (m); or, e.g. (gamma/(gamma-1))*H_s, or C_p*theta_s/g, where C_p = 1000 J/kg. K, etc.  
g_real                   = 9.8
beta_eq                  = 2.2793*(10.**-11.)    		# 2*omega*cos((phi/180)*pi)/R_e, phi=lat, R_e: Earth's radius
U_scale                  = ma.sqrt(g_real*H0_dim)
u_asr                    = 0.45                  		# aspect ration of velocity at sea surface level with respect to mean lower layer
L_d_eq                   = ma.sqrt(U_scale/(beta_eq))	#4262500.        # Equatorial Rossby deformation radius= sqrt(c/beta), where c=sqrt(gH)
R_earth                  = 6371000.+H0_dim       		# Radius of Earth
if topography==1:
   a                     = R_earth/L_d_eq        		# 1.7218 corresponds to H0_dim=10km, #a=Aspect ratio with respect to barotropic equatorial Rossby deformation radius (L_d)
else:
   a                     = 1.7                   		# Earth radius. Aspect ratio with respect to barotropic equatorial Rossby deformation radius, L_d=sqrt(U_bt/(beta))
													    # U_bt=sqrt(g*H0); beta= 2*omega*cos((phi/180)*pi)/R_e; R_e: Earth's radius
T_scale                  = a/2.                  		# Time scale to convert to [day]:e.g.:  2*pi/(beta*L_d_eq)  

delta1                   = 0.6                  		# 0.55#0.57 # H1/H0
delta2                   = 1.0-delta1
p_level1                 = P_s                   		# Pa
p_level2                 = (1.0-delta1)*P_s      		# Pa

# Set the month and year variables
# month_year             = 'January1980'            	# Change this to any desired month and year
# month_year             = 'March1980'             	
month_year               = 'June1980'  
# month_year             = 'September1980'  
# month_year             = 'December1980' 
# ullcr                  = 0.07              			#  (Upper Layer Liquid (water) Content Ratio)

B1                 	     = 1.1563               		# Mean value of theta/theta_s at the lower (first) layer
B2                       = 1.2886               		#

if month_year           == 'March1980':
   t_init_insol          =  60.-spin_up_days  			# days to reach the equilibrium state
   input_folder          = '../Data/Aeolus2.0_Input_ERA5/March_1980' # location of the input files (used for restart of the model) 
   start_file            = 'output_1'            		# output_1.npz
   
elif month_year         == 'June1980':
   t_init_insol          =  152.-spin_up_days			# 30=days to reach the equilibrium state
   input_folder          = '../Data/Aeolus2.0_Input_ERA5/June_1980'
   start_file            = 'output_1'            		# output_1.npz
   
elif month_year         == 'September1980': 
   t_init_insol          =  244.-spin_up_days			# days to reach the equilibrium state
   input_folder          = '../Data/Aeolus2.0_Input_ERA5/September_1980'
   start_file            = 'output_1'            		# output_1.npz
         
elif month_year         == 'December1980':
   t_init_insol          =  335.-spin_up_days			# - days to reach the equilibrium state
   input_folder          = '../Data/Aeolus2.0_Input_ERA5/December_1980'
   start_file            = 'output_1'            		# output_1.npz
   
else:
   t_init_insol          =  0.  
   input_folder          = '../Data/Aeolus2.0_Input_ERA5/January_1980'   # "initial0g"             
   start_file            = 'output_1'            		# output_1.npz
   
B3                       = 1.35                  		# Mean value of g*theta/theta_o at the upper (third) layer, it does not exist in two-layer configuration
#g = g0*np.array([[B1, B1, B1],[B1, B2, B2],[B1,B2,B3]])


if restart_from_output:
	if not os.path.exists(output_folder):
		print('\nfail --> output folder not exist for restart ("{}")'.format(output_folder))
		exit()
	input_folder         = output_folder
	start_file           = sorted(glob.glob('{}/{}*.npz'.format(output_folder, output_file)))[-1].split('/')[-1].replace('.npz','') 

