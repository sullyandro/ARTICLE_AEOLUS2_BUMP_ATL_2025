################################################################################
#
# Author:    Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Coauthors: Masoud Rostami, Stefan Petri
# Date:      04.12.2025
# Type:      Python3 + CDO
#
# Description:
# Script to produce Aeolus 2.0 adimensional input data using Reanalysis.
#
# The main variables used from reanalysis are TA, UA, and VA. 
#
# At the end b1, b2, U1 (u1ph), U2 (u2ph), V1 (u1th), V2 (u2th), are saved in some formats to be read by Aeolus2.
#
# Other variables (h1, h2, q1, q2, ...) are included as a rest state (equal zero).
#
################################################################################

import os
import numpy as np
import math as ma
import warnings
import humanize
from scipy.io import savemat
from NetCDFOutput import NetCDFOutput
from netCDF4 import Dataset, num2date
warnings.filterwarnings('ignore')

def   size(f): return humanize.naturalsize(os.path.getsize(f))
def exists(f): return os.path.exists(f)
    
    
print('')
print('Aeolus2.0 Startup conditions Preparation from Reanalysis')
print('')


# General definitions of dimentional parameters

g                   = 9.8                               # Earth gravity

omega               = 7.2720e-05                        # Earth angular speed (rad/s)

H0_dim              = 10000                             # Scale of (pseudo)-height (m); 10km; https://www.e3s-conferences.org/articles/e3sconf/pdf/2019/02/e3sconf_icst2018_04002.pdf 

R_earth             = 6371000 + H0_dim                  # Earth radius = 6371 km ; plus the adjustment between topograph and atmosphere

U_scale             = ma.sqrt(g*H0_dim)                 # Barotropic velocity speed, 313 m/s, also U scale 

beta_eq             = 2.2793e-11                        # Rossby parameter from Coriolis; 2*omega*cos((phi/180)*pi)/R_earth; phi=lat=0=equator

L_d_eq              = ma.sqrt(U_scale/(beta_eq))        # Equatorial Rossby Deformation Radius; 3706030.236556578 (matlab); 3706003.060829316 (python) 

R_o_Cp              = 0.285                             # Ratio of the universal gas constant (R) to the specific heat capacity at constant pressure (Cp) ; ~0.286 for air

a                   = R_earth/L_d_eq                    # Aspect ratio with respect to barotropic equatorial Rossby deformation radius (L_d) ; 1.7218010604049752

T_scale             = a/2.                              # Time scale to convert from the nondimensinal scale to dimensional [day] ; 0.8609005302024876
                                                        # e.g.: Aeolus2 time (0, 0.28889449, ...) --> Real Earth time (0, 0.24870942, ...) ; the model runs in a "faster time scale"
                                                        # Real Earth time (time) = Aeolus2 time (t) * T_scale  ; the result is a fraction of a day --> [0, 0.24870942, ..., 1 (day), 1.24870942, ...]



################################
#        Aeolus2 grid          # 
################################

# Full path of the grid in netcdf format
# Only the lat and lon are important to be read from this file.

# Aeolus2 nc output files are in format (-180:180, -90:90). This grid sample is in the right way and will be used by cdo ahead as desired final grid.
# Aeolus2 npz files are with diffente setup, explaned ahead. 
                                                        
grid_sample         = '../Data/Aeolus2.0_Grids/aeolus2_grid_768x384.nc'
# latitude  = 384 (-89.64165, -89.17743, -88.71048, ..., 88.71048, 89.17743, 89.64165)
# longitude = 768 (-180.0000, -179.5312, -179.0625, ..., 178.5938, 179.0625, 179.5312)

print()
print('Aeolus2 grid -->', grid_sample, size(grid_sample))





################################
#    Aeolus2 Startup Output    # 
################################

# Aeolus2 Startup Output name (to be add '.nc', '.npz', '.mat')

aeolus2_startup_output     = '../Data/Aeolus2.0_Input_ERA5/Aeolus2_startup_768x384_from_ERA5_6hours_1980-06-01-00_1980-06-30-18_timmean'




################################
#       Reanalysis Data        # 
################################

# Example: ERA-Interim
# latitude  =  241 (90, 89.25, 88.5, ...,  -88.5, -89.25,    -90)
# longitude =  480 ( 0,  0.75,  1.5, ..., 357.75,  358.5, 359.25)

# Example: ERA5
# latitude  =  721 (90, 89.75, 89.5, ...,  -89.5, -89.75,    -90)
# longitude = 1440 ( 0,  0.25, 0.75, ..., 359.25,  359.5, 359.75)

reanalysis          = '../Data/ERA5/ERA5_6hours_1980-06-01-00_1980-06-30-18.nc'

print()
print('Reanalysis   -->', reanalysis, size(reanalysis))

                                                 
# CDO selection parameters

cdo_seldate         = '1980-06-01T00:00:00,1980-06-30T23:00:00'  

cdo_time_operator   = '-timmean'                        # change here to a empty string " " for no time mean application 

cdo_selname         = 't,u,v'

cdo_sellevel        = '250,450,650,800,900,1000'        # hPa

reanalysis_remap    = '../Data/ERA5/{}_remapbil_{}_variables_selection.nc'.format(reanalysis.split('/')[-1][:-3], grid_sample.split('/')[-1][:-3]) # ERA5_6hours_1980-06-01-00_1980-06-30-18_remapbil_aeolus2_grid_768x384_variables_selection.nc


# Data selection and Interpolation of Reanalysis to Aeolus2 grid

cmd                 = 'cdo -P 8 -O remapbil,{0} {1} -sellevel,{2} -selname,{3} -seldate,{4} {5} {6}'.format(grid_sample, cdo_time_operator, cdo_sellevel, cdo_selname, cdo_seldate, reanalysis, reanalysis_remap)
                    
                    # cdo -P 8 -O remapbil,aeolus2_grid_768x384.nc -timmean -sellevel,250,450,650,800,900,1000 -selname,t,u,v -seldate,1980-06-01T00:00:00,1980-06-30T23:00:00 
                    # ERA5_6hours_1980-06-01-00_1980-06-30-18.nc 
                    # ERA5_6hours_1980-06-01-00_1980-06-30-18_remapbil_aeolus2_grid_768x384_variables_selection.nc
print()
print('cmd -->', cmd)

os.system(cmd)

print()
if exists(reanalysis_remap):
    print('done -->', reanalysis_remap, size(reanalysis_remap))



# Reading Regrided Reanalysis 

print()
print('+++++++++++++++++++++++++++++++++')
print('# Reading NetCDF4 file')
print()

if exists(reanalysis_remap): print('input -->', reanalysis_remap, size(reanalysis_remap))
else:                        print('fail --> not found',  reanalysis_remap) ; exit()

dat = Dataset(reanalysis_remap, mode='r')
lon = dat.variables['longitude'][:]
lat = dat.variables['latitude' ][:]
try:    
	lev       = dat.variables['level'][:]
except: 
	lev       = dat.variables['pressure_level'][:]
try:    
	tim       = dat.variables['time'][:]
	tim_units = dat.variables['time'].units
except: 
	tim       = dat.variables['valid_time'][:]
	tim_units = dat.variables['valid_time'].units
dtime = num2date(tim[:], tim_units).astype('datetime64[ms]').astype('O')

variables = list(dat.variables.keys())

tim_size = len(tim)
lat_size = len(lat)
lon_size = len(lon)
lev_size = len(lev)


print()
print('lats: {} [{:2.6f}, {:2.6f}, ..., {:2.6f}, {:2.6f}]'.format(len(lat), lat[0], lat[1], lat[-2], lat[-1])) # 384 [ -89.64165 -89.17743 ... 89.17743 89.64165 ]
print()
print('lons: {} [{:2.6f}, {:2.6f}, ..., {:2.6f}, {:2.6f}]'.format(len(lon), lon[0], lon[1], lon[-2], lon[-1])) # 768 [ -180.0 -179.5312 ... 179.0625 179.5312 ]
print()
print('time:', tim_size, dtime)                                                                                # 1 [datetime.datetime(1980, 1, 1, 9, 0)]
print()
print('level:', lev_size, lev)                                                                                 # 6 [ 250  450  650  800  900 1000]
print()
print(variables)
print()

TA = dat.variables['t'][:] # Temperature         (K)
UA = dat.variables['u'][:] # U component of wind (m s**-1)
VA = dat.variables['v'][:] # V component of wind (m s**-1)

print()
print('TA -->', TA.shape, dat.variables['t'].dimensions) # (1, 7, 384, 768) ('time', 'level', 'latitude', 'longitude')
print('UA -->', TA.shape, dat.variables['u'].dimensions) # (1, 7, 384, 768) ('time', 'level', 'latitude', 'longitude')
print('VA -->', TA.shape, dat.variables['v'].dimensions) # (1, 7, 384, 768) ('time', 'level', 'latitude', 'longitude')
print()

dat.close()
os.system('rm -rf '+ reanalysis_remap)



####################################
# Reanalysis Conversion to Aeolus2 # 
####################################


if np.sort(np.array(lev,dtype=int)).tolist() == [250, 450, 650, 800, 900, 1000]:	
	
	print('')
	print('LEVELS CASE -->', [250, 450, 650, 800, 900, 1000])
	
	# Variables by level

	TA_250  = TA[:,lev==250, :,:]
	TA_450  = TA[:,lev==450, :,:]
	TA_650  = TA[:,lev==650, :,:]
	TA_800  = TA[:,lev==800, :,:]
	TA_900  = TA[:,lev==900, :,:]
	TA_1000 = TA[:,lev==1000,:,:]

	UA_250  = UA[:,lev==250, :,:]
	UA_450  = UA[:,lev==450, :,:]
	UA_650  = UA[:,lev==650, :,:]
	UA_800  = UA[:,lev==800, :,:]
	UA_900  = UA[:,lev==900, :,:]
	UA_1000 = UA[:,lev==1000,:,:]

	VA_250  = VA[:,lev==250, :,:]
	VA_450  = VA[:,lev==450, :,:]
	VA_650  = VA[:,lev==650, :,:]
	VA_800  = VA[:,lev==800, :,:]
	VA_900  = VA[:,lev==900, :,:]
	VA_1000 = VA[:,lev==1000,:,:]


	# Convert Temperature to Potential Temperature (PotTA)

	PotTA_250  = TA_250  * ( (1000/250 )**R_o_Cp )
	PotTA_450  = TA_450  * ( (1000/450 )**R_o_Cp )
	PotTA_650  = TA_650  * ( (1000/650 )**R_o_Cp )
	PotTA_800  = TA_800  * ( (1000/800 )**R_o_Cp )
	PotTA_900  = TA_900  * ( (1000/900 )**R_o_Cp )
	PotTA_1000 = TA_1000 * ( (1000/1000)**R_o_Cp )


	# Compute variables in nondimensional scale
	
	# Masoud note:
	# If the pressure levels are not equi-distanced, then the distance needs to be included
	
	UA_upp      = ( (UA_250*2.5 + UA_450*2.0)/4.5 ) / U_scale                                 

	UA_low      = ( (UA_650*2.0 + UA_800*1.5 + UA_900*1.0 + UA_1000*1.0)/5.5 ) / U_scale


	VA_upp      = ( (VA_250*2.5 + VA_450*2.0)/4.5 ) / U_scale                                 

	VA_low      = ( (VA_650*2.0 + VA_800*1.5 + VA_900*1.0 + VA_1000*1.0)/5.5 ) / U_scale


	PotTA_upp   = ( (PotTA_250*2.5 + PotTA_450*2.0)/4.5 )                             

	PotTA_low   = ( (PotTA_650*2.0 + PotTA_800*1.5 + PotTA_900*1.0 + PotTA_1000*1.0)/5.5 )

else:
		
	print('')
	print('LEVELS CASE --> No calculations prepared for the given levels! Code needs to be updated.')
	print('')
	exit()
	

# Buoyancy calculation in nondimensional scale

# Masoud notes:
# To compare comparables B2, TS, should be the same for all months
# Note that b(dimensional) = g*theta/theta_s=g*b', b' is non-dim.
# TS = min(PotTA_low) ; 1st_jan=250.0277 ; monthly_aver=262.1009

TS          = np.min(PotTA_low)     # 244.73145 for era5 1980-06-01-00_1980-06-30-18_timmean

print()
print('TS -->', TS)
print()
	
B1_mean     = np.mean(PotTA_low) / TS
B2_mean     = np.mean(PotTA_upp) / TS

B_low       = ( PotTA_low / TS ) - B1_mean
B_upp       = ( PotTA_upp / TS ) - B2_mean



# Renaming or defining variables to the final output

info = {}

# Essentials

b1          = B_low         ;  info['b1'       ] = {'data':b1       , 'longname':'Buoyancy layer 1'                            , 'units':'[g*theta/theta_s]'                       }
b2          = B_upp         ;  info['b2'       ] = {'data':b2       , 'longname':'Buoyancy layer 2'                            , 'units':'[g*theta/theta_s]'                       }

u1ph        = UA_low        ;  info['u1ph'     ] = {'data':u1ph     , 'longname':'Zonal (azimuthal) velocity layer 1'          , 'units':'[(g*H)^0.5], [beta*L_d^2] at the equator'}
u2ph        = UA_upp        ;  info['u2ph'     ] = {'data':u2ph     , 'longname':'Zonal (azimuthal) velocity layer 2'          , 'units':'[(g*H)^0.5], [beta*L_d^2] at the equator'}

u1th        = VA_low        ;  info['u1th'     ] = {'data':u1th     , 'longname':'Meridional velocity layer 1'                 , 'units':'[(g*H)^0.5], [beta*L_d^2] at the equator'}
u2th        = VA_upp        ;  info['u2th'     ] = {'data':u2th     , 'longname':'Meridional velocity layer 2'                 , 'units':'[(g*H)^0.5], [beta*L_d^2] at the equator'}

# Resting state (zero)

zeros       = b1*0

h1          = zeros         ;  info['h1'       ] = {'data':h1       , 'longname':'Pseudo-height of layer 1'                    , 'units':'[H]'                                     }
h2          = zeros         ;  info['h2'       ] = {'data':h2       , 'longname':'Pseudo-height of layer 2'                    , 'units':'[H]'                                     }

q1          = zeros         ;  info['q1'       ] = {'data':q1       , 'longname':'Bulk of Specific humidity at layer 1'        , 'units':'[(L_v.g/(C_p.theta_s))Kg/Kg]'            }
q2          = zeros         ;  info['q2'       ] = {'data':q2       , 'longname':'Bulk of Specific humidity at layer 2'        , 'units':'[(L_v.g/(C_p.theta_s))Kg/Kg]'            }
w1          = zeros         ;  info['w1'       ] = {'data':w1       , 'longname':'Bulk of Precipitable Water at layer 1'       , 'units':'[(L_v.g/(C_p.theta_s))Kg/Kg]'            }
w2          = zeros         ;  info['w2'       ] = {'data':w2       , 'longname':'Bulk of Precipitable Water at layer 1'       , 'units':'[(L_v.g/(C_p.theta_s))Kg/Kg]'            }
Prec1       = zeros         ;  info['Prec1'    ] = {'data':Prec1    , 'longname':'Bulk of Precipitaion at layer 1'             , 'units':'[(L_v.g/(C_p.theta_s))Kg/Kg]'            }

CC1         = zeros         ;  info['CC1'      ] = {'data':CC1      , 'longname':'CLWC (Cloud Liquid Water Content) at layer 1', 'units':'[Q/T]'                                   }
DD1         = zeros         ;  info['DD1'      ] = {'data':DD1      , 'longname':'Downdraft to layer 1 (balanced)'             , 'units':'[Q/T]'                                   }
DDr         = zeros         ;  info['DDr'      ] = {'data':DDr      , 'longname':'Downdraft to layer 1 (unbalanced)'           , 'units':'[Q/T]'                                   }
Ev          = zeros         ;  info['Ev'       ] = {'data':Ev       , 'longname':'Sea surface evaporation'                     , 'units':'[Q/T]'                                   }

RT1         = zeros         ;  info['RT1'      ] = {'data':RT1      , 'longname':'Radiative transfer flux at layer 1'          , 'units':'[HB]'                                    }
RT2         = zeros         ;  info['RT2'      ] = {'data':RT2      , 'longname':'Radiative transfer flux at layer 2'          , 'units':'[HB]'                                    }

menthalpy   = zeros         ;  info['menthalpy'] = {'data':menthalpy, 'longname':'Menthalpy = h1+H1-ep1*q1-ep1*Q01'            , 'units':'[unknown units]'                         }

 

# Checking

print()
print('Checking variables output before saving')
print()
print( 'time      -> ', tim.shape)
print( 'latitude  -> ', lat.shape)
print( 'longitude -> ', lon.shape)
print()
print( 'b1        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   b1.shape, np.min(b1   ), np.mean(b1   ), np.max(b1   )) )  # min =  0.0000  |  mean =  0.1638  |  max = 0.2707
print( 'b2        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   b2.shape, np.min(b2   ), np.mean(b2   ), np.max(b2   )) )  # min =  0.0480  |  mean =  0.1638  |  max = 0.2591
print( 'u1ph      ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format( u1ph.shape, np.min(u1ph ), np.mean(u1ph ), np.max(u1ph )) )  # min = -0.0702  |  mean =  0.0074  |  max = 0.0896
print( 'u2ph      ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format( u2ph.shape, np.min(u2ph ), np.mean(u2ph ), np.max(u2ph )) )  # min = -0.0962  |  mean =  0.0345  |  max = 0.2171
print( 'u1th      ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format( u1th.shape, np.min(u1th ), np.mean(u1th ), np.max(u1th )) )  # min = -0.0654  |  mean = -0.0002  |  max = 0.0692
print( 'u2th      ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format( u2th.shape, np.min(u2th ), np.mean(u2th ), np.max(u2th )) )  # min = -0.1145  |  mean =  0.0002  |  max = 0.1411
print()
print( 'h1        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   h1.shape, np.min(h1   ), np.mean(h1   ), np.max(h1   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'h2        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   h2.shape, np.min(h2   ), np.mean(h2   ), np.max(h2   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'q1        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   q1.shape, np.min(q1   ), np.mean(q1   ), np.max(q1   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'q2        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   q2.shape, np.min(q2   ), np.mean(q2   ), np.max(q2   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'w1        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   w1.shape, np.min(w1   ), np.mean(w1   ), np.max(w1   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'w2        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   w2.shape, np.min(w2   ), np.mean(w2   ), np.max(w2   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'Prec1     ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(Prec1.shape, np.min(Prec1), np.mean(Prec1), np.max(Prec1)) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'CC1       ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(  CC1.shape, np.min(CC1  ), np.mean(CC1  ), np.max(CC1  )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'DD1       ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(  DD1.shape, np.min(DD1  ), np.mean(DD1  ), np.max(DD1  )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'DDr       ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(  DDr.shape, np.min(DDr  ), np.mean(DDr  ), np.max(DDr  )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'Ev        ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(   Ev.shape, np.min(Ev   ), np.mean(Ev   ), np.max(Ev   )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'RT1       ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(  RT1.shape, np.min(RT1  ), np.mean(RT1  ), np.max(RT1  )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'RT2       ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(  RT2.shape, np.min(RT2  ), np.mean(RT2  ), np.max(RT2  )) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000
print( 'menthalpy ->  {}  |  min = {:2.4f}  |  mean = {:2.4f}  |  max = {:2.4f}'.format(menthalpy.shape, np.min(menthalpy), np.mean(menthalpy), np.max(menthalpy)) )  # min =  0.0000  |  mean =  0.0000  |  max = 0.0000





# Saving as NetCDF4

save_nc = 1

if save_nc == 1: 

    print()
    print()
    print('+++++++++++++++++++++++++++++++++')
    print('# Saving NetCDF4 file')
    print()
    
    fout_nc = aeolus2_startup_output + '.nc'

    foo = Dataset(fout_nc, 'w', format='NETCDF4_CLASSIC')

    foo.createDimension('time', None)
    foo.createDimension('latitude',  len(lat))
    foo.createDimension('longitude', len(lon))

    lats                = foo.createVariable('latitude', 'f4', ('latitude'), zlib=True)
    lats.units          = 'degrees_north'
    lats.long_name      = 'latitude'
    lats.axis           = 'Y'
    lats[:]             = lat[:]

    lons                = foo.createVariable('longitude', 'f4', ('longitude'), zlib=True)
    lons.units          = 'degrees_east'
    lons.long_name      = 'longitude'
    lons.axis           = 'X'
    lons[:]             = lon[:]

    times               = foo.createVariable('time', 'f4', ('time'), zlib=True)
    times.units         = tim_units
    times.calendar      = 'standard'
    times.standard_name = 'time'
    times[:]            = tim[-1]

    for ncvar in info:  # ['b1', 'b2', ...]
        
        globals()['ncvar_' + ncvar]               = foo.createVariable(ncvar, float, ('time', 'latitude', 'longitude'), zlib=True)
        globals()['ncvar_' + ncvar].units         = info[ncvar]['units']
        globals()['ncvar_' + ncvar].long_name     = info[ncvar]['longname']
        globals()['ncvar_' + ncvar].missing_value = np.nan
        globals()['ncvar_' + ncvar][:]            = info[ncvar]['data'][:]
        
    foo.comment         = 'Aeolus 2.0 startup data (nondimensional scale) from Reanalysis ({})'.format(reanalysis.split('/')[-1])
    foo.close()

    if exists(fout_nc):
        print('done -->', fout_nc, size(fout_nc))




# Saving as npz

save_npz = 1

if save_npz == 1: 
    
    print()
    print()
    print('+++++++++++++++++++++++++++++++++')
    print('# Saving NPZ file')
    print()
    
    fout_npz = aeolus2_startup_output + '.npz'

    # Creating t, lamda, and theta to match with npz samples

    """
    Given the goal here be to prepare a single timestep for Aeolus initialization, the first step will be considered as 0. Thus, t = 0.0 / T_scale 

    For a longer time axis conversion (not the actual case), see below.
    e.g.: Era5 6hs frequency data --> [0h, 6h, 12h, 18h] --> in day fractions --> [0., 0.333, 0.666, 1.] --> divide by "Aeolus2 T_scale 0.8609005302024876" --> t = [0., 0.38719146, 0.77438292, 1.16157438]
    """
    
    # Time here has to be in day fractions
    
    time  = np.array([0.0])         				# time in Real Earth time scale
    t     = time / T_scale          				# time in nondimensional scale for Aeolus2


    # Latitude here has to be 90:-90

    if lat[0]<0: # shifting from -90:90 to 90:-90
        
        theta = (lat[::-1] + 90) / (180/np.pi)      # latitude in nondimensional scale for Aeolus2

        for npzvar in info:  # ['b1', 'b2', ...]
        
            globals()[npzvar] = globals()[npzvar][:, ::-1, :]
    
    else:
        
        theta = (lat[:] + 90) / (180/np.pi)         # latitude in nondimensional scale for Aeolus2
    
    
    # longitude here has to be 0:360
    
    if np.any(lon<0): # shifting from -180:180 to 0:360
            
        lon_shift = np.concatenate([np.where(lon>=0)[0], np.where(lon<0)[0]], axis=0)
                
        lon_0_360 = np.where(lon<0, lon+360, lon)[lon_shift]
        
        lamda = (lon_0_360) / (180/np.pi)           # longitude in nondimensional scale for Aeolus2
        
        
        # For some reason the data matrix in npz files are not shifted to 0:360. 
        # But, below is the way to do the shifting, in case needed in the future
        
        # for npzvar in info:  # ['b1', 'b2', ...]
        
            # globals()[npzvar] = globals()[npzvar][:, :, lon_shift]
    
    else:
        
        lamda = (lon) / (180/np.pi)                 # longitude in nondimensional scale for Aeolus2
    
    
    print('theta and lamda are the latitude and longitude in nondimensional scale for Aeolus2')
    print()
    print('theta (lat 90:-90): {} [{:2.6f}, {:2.6f}, ..., {:2.6f}, {:2.6f}]'.format(len(theta), theta[0], theta[1], theta[-2], theta[-1])) # 384 [3.135338, 3.127236, ..., 0.014357, 0.006254]
    print()
    print('lamda (lon  0:360): {} [{:2.6f}, {:2.6f}, ..., {:2.6f}, {:2.6f}]'.format(len(lamda), lamda[0], lamda[1], lamda[-2], lamda[-1])) # 768 [0.000000, 0.008181, ..., 6.266823, 6.275004]
    print()


    # Aeolus2 Npz files contains only 1 timestep, and in (lon, lat) shape, so removing time dimension and swaping axis 
    
    for npzvar in info:  # ['b1', 'b2', ...]
    
        globals()[npzvar] = np.swapaxes(np.squeeze(globals()[npzvar][:, :, :]), 0, 1)  # from [lat,lon] to [lon,lat]
        
        
    # Saving
    np.savez_compressed(fout_npz, t=t, theta=theta, lamda=lamda, 
                        b1=b1, b2=b2, u1ph=u1ph, u2ph=u2ph, u1th=u1th, u2th=u2th, h1=h1, h2=h2, q1=q1, q2=q2, w1=w1, w2=w2, 
                        Prec1=Prec1, CC1=CC1, DD1=DD1, DDr=DDr, Ev=Ev, RT1=RT1, RT2=RT2, menthalpy=menthalpy)

    if exists(fout_npz):
        print('done -->', fout_npz, size(fout_npz))




# Saving as mat

save_mat = 1

if save_mat == 1: 
    
    print()
    print()
    print('+++++++++++++++++++++++++++++++++')
    print('# Saving mat file')
    print()
    
    fout_mat = aeolus2_startup_output + '.mat'
        
    if exists(fout_npz): 
        print('input -->', fout_npz, size(fout_npz))
    else:                        
        print('fail --> not found', fout_npz) ; exit()

    npz_loaded = np.load(fout_npz)
    
    savemat(fout_mat, npz_loaded)

    if exists(fout_mat):
        print('done  -->', fout_mat, size(fout_mat))
        
        
        



# End

