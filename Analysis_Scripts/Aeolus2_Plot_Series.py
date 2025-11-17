################################################################################
#
# Author:       Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Colaborators: Masoud Rostami
# Date:         01.08.2025
# Type:         Python3
#
# Description:
# Script to plot graphics for Aeolus2 analysis of Bump evolution.
#
################################################################################


import os
import sys
import glob
import time
import argparse
import humanize
import warnings
import numpy as np
import matplotlib as mpl #; mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import cartopy.crs as ccrs
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from scipy.ndimage import uniform_filter
plt.style.use('classic')
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('--remap', help='data from remapbil', action='store_true')
args = parser.parse_args()

remap = args.remap

def   size(f): return humanize.naturalsize(os.path.getsize(f))
def exists(f): return os.path.exists(f)
    

print('')
print('')
print('Aeolus2.0 Analysis Preparation - Plot Series')


cases = {
'Dry_Barotropic_Weak'  :{'typ':'Barotropic', 'bump':'weak',   'dir':'../Data/Aeolus2.0_Output_Dry_Barotropic_Weak',   'nc':'Aeolus2.0_Output_Dry_Barotropic_Weak_minus_Dry_Control.nc'  },
'Dry_Barotropic_Strong':{'typ':'Barotropic', 'bump':'strong', 'dir':'../Data/Aeolus2.0_Output_Dry_Barotropic_Strong', 'nc':'Aeolus2.0_Output_Dry_Barotropic_Strong_minus_Dry_Control.nc'},
'Dry_Baroclinic_Weak'  :{'typ':'Baroclinic', 'bump':'weak',   'dir':'../Data/Aeolus2.0_Output_Dry_Baroclinic_Weak',   'nc':'Aeolus2.0_Output_Dry_Baroclinic_Weak_minus_Dry_Control.nc'  },
'Dry_Baroclinic_Strong':{'typ':'Baroclinic', 'bump':'strong', 'dir':'../Data/Aeolus2.0_Output_Dry_Baroclinic_Strong', 'nc':'Aeolus2.0_Output_Dry_Baroclinic_Strong_minus_Dry_Control.nc'},
'MC_Barotropic_Weak'   :{'typ':'Barotropic', 'bump':'weak',   'dir':'../Data/Aeolus2.0_Output_MC_Barotropic_Weak',    'nc':'Aeolus2.0_Output_MC_Barotropic_Weak_minus_MC_Control.nc'    },
'MC_Barotropic_Strong' :{'typ':'Barotropic', 'bump':'strong', 'dir':'../Data/Aeolus2.0_Output_MC_Barotropic_Strong',  'nc':'Aeolus2.0_Output_MC_Barotropic_Strong_minus_MC_Control.nc'  },
'MC_Baroclinic_Weak'   :{'typ':'Baroclinic', 'bump':'weak',   'dir':'../Data/Aeolus2.0_Output_MC_Baroclinic_Weak',    'nc':'Aeolus2.0_Output_MC_Baroclinic_Weak_minus_MC_Control.nc'    },
'MC_Baroclinic_Strong' :{'typ':'Baroclinic', 'bump':'strong', 'dir':'../Data/Aeolus2.0_Output_MC_Baroclinic_Strong',  'nc':'Aeolus2.0_Output_MC_Baroclinic_Strong_minus_MC_Control.nc'  },
}

reg = 'Atlantic'
loc = [-55, -35, 10, 30]

lat_0 =  20.
lon_0 = -40.

data = { }

for ie, exp in enumerate(cases):

    print()
    print()
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print()
    print('Running -->', exp)
    print()
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print()
    
    
    exp_nc   = cases[exp]['nc'  ]
    exp_dir  = cases[exp]['dir' ]
    exp_typ  = cases[exp]['typ' ]
    exp_bump = cases[exp]['bump']

    data[exp] = {}

    for var in ['u1ph', 'u1th', 'u2ph', 'u2th', 'h1', 'h2', 'b1', 'b2', 'CC1', 'w2']:

        if remap == 0: fin = '{}/{}'.format(exp_dir, exp_nc).replace('.nc','_{}.nc'.format(var))
        if remap == 1: fin = '{}/{}'.format(exp_dir, exp_nc).replace('.nc','_{}_remapbil_1dg.nc'.format(var))

        if exists(fin): print('input -->', fin, size(fin))
        else:           print('fail --> not found',  fin) ; exit()

        print()
        print('+++++++++++++++++++++++++++++++++')
        print('# Reading NetCDF4 file')
        
        dat = Dataset(fin, mode='r')
        lon = dat.variables['longitude'][:]
        lat = dat.variables['latitude' ][:]
        tim = dat.variables['time'][:]
        tim_units = dat.variables['time'].units

        variables = list(dat.variables.keys())
        variables.remove('longitude')
        variables.remove('latitude')
        variables.remove('time')

        tim_size = len(tim)
        lat_size = len(lat)
        lon_size = len(lon)

        # Time variables
         
        dtime = num2date(tim[:], tim_units)
        dtime_hov = dtime.astype('datetime64[ms]').astype('O') # for Hovmoeller plot
        # print(str(dtime[21])) # bump starts at step 21

        # Creating a time var for ploting
        times = []
        tim_range = range(tim_size)
        for i in tim_range:
            times.append( str(dtime[i])[:10] )
            # ['1980-06-01', '1980-06-01', ..., '1980-09-10', '1980-09-11']


        print()
        print(variables)
        print()
        print('lats:', lat_size, '[', lat[0], lat[1], '...', lat[-2], lat[-1], ']')
        print()
        print('lons:', lon_size, '[', lon[0], lon[1], '...', lon[-2], lon[-1], ']')
        print()
        print('time:', tim_size, '[', dtime[0], dtime[1], '...', dtime[-2], dtime[-1], ']')


        # Data variables 
        
        if var == 'u1ph': data[exp]['u1_t' ] = dat.variables['u1ph'][:,:,:] # Zonal (azimuthal) velocity layer 1           ; units: [(g*H)^0.5], [beta*L_d^2] at the equator
        if var == 'u1th': data[exp]['v1_t' ] = dat.variables['u1th'][:,:,:] # Meridional velocity layer 1                  ; units: [(g*H)^0.5], [beta*L_d^2] at the equator
        
        if var == 'u2ph': data[exp]['u2_t' ] = dat.variables['u2ph'][:,:,:] # Zonal (azimuthal) velocity layer 2           ; units: [(g*H)^0.5], [beta*L_d^2] at the equator
        if var == 'u2th': data[exp]['v2_t' ] = dat.variables['u2th'][:,:,:] # Meridional velocity layer 2                  ; units: [(g*H)^0.5], [beta*L_d^2] at the equator

        if var == 'h1':   data[exp]['h1_t' ] = dat.variables['h1'  ][:,:,:] # Pseudo-height layer 1                        ; units: H
        if var == 'h2':   data[exp]['h2_t' ] = dat.variables['h2'  ][:,:,:] # Pseudo-height layer 2                        ; units: H
        
        if var == 'b1':   data[exp]['b1_t' ] = dat.variables['b1'  ][:,:,:] # Buoyancy layer 1                             ; units: g*theta/theta_s
        if var == 'b2':   data[exp]['b2_t' ] = dat.variables['b2'  ][:,:,:] # Buoyancy layer 2                             ; units: g*theta/theta_s
        
        if var == 'CC1':  data[exp]['cc1_t'] = dat.variables['CC1' ][:,:,:] # CLWC (Cloud Liquid Water Content) at layer 1 ; units: Q/T
        if var == 'w2':   data[exp]['w2_t' ] = dat.variables['w2'  ][:,:,:] # Bulk of Precipitable Water at layer 2        ; units: (L_v.g/(C_p.theta_s))Kg/Kg

        dat.close()
    
            
print()
print('+++++++++++++++++++++++++++++++++')
print('# General constants and definitions')
print()


grilly, grillx      = np.meshgrid(lat, lon)

Nx                  = lon_size  
Ny                  = lat_size      

theta_rad           = np.pi/2 - (lat/360) * (2*np.pi)

lamda_rad           = ( (lon+180)*(2*np.pi) )/360


diff_theta_rad      = np.diff(theta_rad)        # u = u_r, v = v_phi
diff_lamda_rad      = np.diff(lamda_rad[::-1])          


##########################################
##########################################

g                   = 9.8                        
H0                  = 10000                         # https://www.e3s-conferences.org/articles/e3sconf/pdf/2019/02/e3sconf_icst2018_04002.pdf
U_bt                = np.sqrt(g*H0)                 # 313ms^-1, also U scale 
omega               = 7.2720e-05                    # Earth angular speed (rad/s)
phi                 = 0                             # degree lat
R_earth             = 6371000+10000                 # 6,371 km, Earth radius (m)
beta                = 2*omega*np.cos((phi/180)*np.pi)/R_earth           # 1.3738e-11 at 53deg
beta40              = 2*omega*np.cos(( 40/180)*np.pi)/R_earth     
L_d_eq              = np.sqrt(U_bt/(beta))                              # 4228 km  # 3380 km
L_d                 = np.sqrt(U_bt/(2*omega*np.sin((phi/180)*np.pi)))   # 1342 km 
LL                  = np.sqrt(U_bt/(2*omega                     ))       
Lambs_p             = (2*omega*R_earth)**2/(g*H0)
U_scale             = U_bt                      
R_gas               = 8.3145                        # J/(K.mol)
R_e_non             = R_earth/L_d_eq            
energy_scale        = (beta**3)*L_d_eq**5       
xexp                = grillx[:,1]               
yexp                = grilly[1,:]               

Ld_scaling          = 1 
Lamb_scaling        = 0 

if Ld_scaling == 1:

    L_scale         = L_d_eq  
    T_scale         = 2*np.pi*R_e_non/(2*omega)     # 2*pi/(beta*L_d_eq) ; s/rad, is equal to 2*pi*R_e_non/(2*omega)
    T_scale_day     = T_scale/86400              
    r_sphere        = R_earth/L_d_eq            
    degree_lon      = 2*np.pi*r_sphere*L_d_eq/360   
    T_scale_d       = 2*np.pi/(beta*L_d_eq)/86400   
    
if Lamb_scaling == 1:

    T_scale         = 2*np.pi/(omega)               
    T_scale_day     = T_scale/86400             
    L_scale         = R_earth                   
    H_coeff         = 1/np.sqrt(Lambs_p)            
    r_sphere        = R_earth/L_scale           
    
    
Tfactor             = T_scale_day               
surf                = 4*np.pi*r_sphere**2           

# ---- dimensional bulk of specific humidity over H -----

H_q                 = 10000                      
lat_h               = 2.5  *10**6                   # J/Kg
C_p                 = 4.218*10**3                   # valid for water not atm., J/Kg.C
q_real              = 20*0.001                      # (Kg/Kg)
thet_s              = 280                           # K
bulk_q_dim          = lat_h*g*H_q/(C_p*thet_s)      
coeff_q_o_h         = lat_h*g    /(C_p*thet_s)      

# fixed parameters
H1                  = 0.35                       
H2                  = 1-H1                       
B1                  = 1                         
B2                  = 1.1                       


T_s                 = 280                           # K
Lr                  = 9.8                           # degree/Km
h_s                 = np.arange(0,6501,100)                 
Tatm                = (T_s - Lr * h_s/1000)         
f1                  = 1 - (  (np.exp(-g*h_s/(R_gas*Tatm)))  )**0.28 

RT_scale2           = g*beta*(L_d_eq)             /T_s      

RT_scale            = g*H0  *(L_d_eq**3 * beta**3)/T_s 

##########################################
##########################################




print()
print('+++++++++++++++++++++++++++++++++')
print('# Plots Series production')
print()




if remap == 0: dir_figs = '../Figures/Plot_Series'
if remap == 1: dir_figs = '../Figures_remapbil_1dg/Plot_Series'

if not exists(dir_figs): os.makedirs(dir_figs)




    
plots_div_layer_1 = 1


if plots_div_layer_1 == 1:
    
    
    print('')
    print('> Divergence in Layer 1')
    print()


    for ie, exp in enumerate(cases):

        print('exp -->', exp)
        
        ##################################################################
        # Calculation of Divergence
 
        u1 = data[exp]['u1_t' ] # u1_t[:,:,:]   # u1ph
        v1 = data[exp]['v1_t' ] # v1_t[:,:,:]   # u1th


        #---------------------------------
        
        v_sin_thet = np.zeros((tim_size, lat_size, lon_size))

        for nn in range(0, lat_size):  # theta

            v_sin_thet[:, nn, :] = v1[:, nn, :] * np.sin(theta_rad[nn])  
            
        #---------------------------------
             
        div1 = np.zeros((tim_size, lat_size-1, lon_size-1))

        for ii in range(0, lon_size-1):
            
            div1[:, :, ii] = np.diff(v_sin_thet[:, :, ii]) / diff_theta_rad[:]
        
        #---------------------------------

        div2 = np.zeros((tim_size, lat_size-1, lon_size-1))

        for jj in range(0, lat_size-1):
            
            div2[:, jj, :] = (1/np.sin(theta_rad[jj])) * (np.diff(u1[:, jj, :]) / diff_lamda_rad[::-1])  # 0<theta<pi

        #---------------------------------

        div3 = np.zeros((tim_size, lat_size-1, lon_size-1))
        
        for jj in range(0, lat_size-1):
            
            div3[:, jj, :] = (1/np.sin(theta_rad[jj])) * (div1[:, jj, :])   # 0<theta<pi

        #---------------------------------
            
        div4            = (1/r_sphere) * ( div2 + div3 ) 
        
        div4[:,  0,  :] = 0.5*( div4[:, -2,  :] + div4[:, 1, :] ) 
        
        div4[:,  :,  0] = 0.5*( div4[:,  :, -2] + div4[:, :, 1] )
        
        div4[:, -1,  :] = 0.5*( div4[:, -2,  :] + div4[:, 1, :] ) 
        
        div4[:,  :, -1] = 0.5*( div4[:,  :, -2] + div4[:, :, 1] )

        #--------------------------------- 
        
        data[exp]['div_layer_1_sum_glb'] = np.nansum(np.abs(div4[:, :, :]), axis=(1,2))
        
        data[exp]['div_layer_1_plt_glb'] = (data[exp]['div_layer_1_sum_glb'])/1000	                # At the legend is added "(x 10^3)"
        
        #--------------------------------- 
    
        # Meridional mean 
        
        data[exp]['div_layer_1_mean_glb_inlat'] = np.nanmean(np.abs(div4[:, :, :]), axis=(0,1))*100	# At the legend is added "(x 10^-2)"
        
        # Zonal mean
        
        data[exp]['div_layer_1_mean_glb_inlon'] = np.nanmean(np.abs(div4[:, :, :]), axis=(0,2))*100	# At the legend is added "(x 10^-2)"
        
        #--------------------------------- 
        
        ##################################################################
        
    
    
    print()
    print('# Ploting -->  Sum of Divergence  =  Wave activity  (Layer 1)  -->  Simple line plot')
    print()
            
    fig, ax = plt.subplots(figsize=(10,5))

    plt.grid(True)
    
    plt.plot(tim_range, data['Dry_Baroclinic_Weak'  ]['div_layer_1_plt_glb'], '--', color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(tim_range, data['Dry_Baroclinic_Strong']['div_layer_1_plt_glb'], '--', color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Weak'  ]['div_layer_1_plt_glb'], '--', color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Strong']['div_layer_1_plt_glb'], '--', color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['div_layer_1_plt_glb'],       color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['div_layer_1_plt_glb'],       color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['div_layer_1_plt_glb'],       color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['div_layer_1_plt_glb'],       color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=12)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    ax.set_ylim(-0.1, 5.5)
    ax.yaxis.set_tick_params(labelsize=12)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.title(' Aeolus2.0 - Wave Activity in Layer 1 ', fontsize=15, pad=15)
    plt.ylabel('(x $10^{3}$)')

    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Wave_Activity_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
        
    
        
    print()
    print('# Ploting -->  Meridional Mean of Divergence  (Layer 1)  -->  Simple line plot')
    print()
                
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)

    lon_plt = lon[:len(data['Dry_Baroclinic_Weak']['div_layer_1_mean_glb_inlat'])]

    plt.plot(lon_plt, data['Dry_Baroclinic_Weak'  ]['div_layer_1_mean_glb_inlat'], color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Baroclinic_Strong']['div_layer_1_mean_glb_inlat'], color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Weak'  ]['div_layer_1_mean_glb_inlat'], color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Strong']['div_layer_1_mean_glb_inlat'], color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Weak'   ]['div_layer_1_mean_glb_inlat'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Strong' ]['div_layer_1_mean_glb_inlat'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Weak'   ]['div_layer_1_mean_glb_inlat'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Strong' ]['div_layer_1_mean_glb_inlat'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    plt.xlabel('Longitude')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_xticks([-180, -135,  -90,  -45,    0,   45,   90,  135,  180])
    ax.set_xticklabels([u'180\N{DEGREE SIGN}', u'135\N{DEGREE SIGN}W', u'90\N{DEGREE SIGN}W', u'45\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'45\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E', u'135\N{DEGREE SIGN}E', u'180\N{DEGREE SIGN}'])
    ax.set_xlim(-180,180.05)
    
    xx = loc[0:2] # [-55, -35]
    y1, y2 = ax.get_ylim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.title(' Aeolus2.0 - Meridional mean of Divergence in Layer 1 ', fontsize=15, pad=15)
    plt.ylabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Divergence_Meridional_Mean_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    

    
    print()
    print('# Ploting -->  Zonal Mean of Divergence  (Layer 1)  -->  Simple line plot')
    print()
            
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)

    lat_plt = lat[:len(data['Dry_Baroclinic_Weak']['div_layer_1_mean_glb_inlon'])]

    plt.plot(data['Dry_Baroclinic_Weak'  ]['div_layer_1_mean_glb_inlon'], lat_plt, color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(data['Dry_Baroclinic_Strong']['div_layer_1_mean_glb_inlon'], lat_plt, color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(data['Dry_Barotropic_Weak'  ]['div_layer_1_mean_glb_inlon'], lat_plt, color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(data['Dry_Barotropic_Strong']['div_layer_1_mean_glb_inlon'], lat_plt, color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(data['MC_Baroclinic_Weak'   ]['div_layer_1_mean_glb_inlon'], lat_plt, color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(data['MC_Baroclinic_Strong' ]['div_layer_1_mean_glb_inlon'], lat_plt, color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(data['MC_Barotropic_Weak'   ]['div_layer_1_mean_glb_inlon'], lat_plt, color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(data['MC_Barotropic_Strong' ]['div_layer_1_mean_glb_inlon'], lat_plt, color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    plt.ylabel('Latitude')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
    ax.set_ylim(-90,90.05)
    
    y1 = [loc[2], loc[2]] # [10, 10]
    y2 = [loc[3], loc[3]] # [30, 30]
    xx = ax.get_xlim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')  
    
    plt.title(' Aeolus2.0 - Zonal mean of Divergence in Layer 1 ', fontsize=15, pad=15)
    plt.xlabel('(x $10^{-2}$)')

    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Divergence_Zonal_Mean_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)

    
    # End of plots_div_layer_1
    


    
plots_div_layer_2 = 1


if plots_div_layer_2 == 1:
    
    
    print('')
    print('> Divergence in Layer 2')
    print()


    for ie, exp in enumerate(cases):

        print('exp -->', exp)
        
        ##################################################################
        # Calculation of Divergence

        u2 = data[exp]['u2_t' ] # u2_t[:,:,:]   # u2ph
        v2 = data[exp]['v2_t' ] # v2_t[:,:,:]   # u2th


        #---------------------------------
        
        v_sin_thet = np.zeros((tim_size, lat_size, lon_size))

        for nn in range(0, lat_size):  # theta

            v_sin_thet[:, nn, :] = v2[:, nn, :] * np.sin(theta_rad[nn])  
            
        #---------------------------------
             
        div1 = np.zeros((tim_size, lat_size-1, lon_size-1))

        for ii in range(0, lon_size-1):
            
            div1[:, :, ii] = np.diff(v_sin_thet[:, :, ii]) / diff_theta_rad[:]
        
        #---------------------------------

        div2 = np.zeros((tim_size, lat_size-1, lon_size-1))

        for jj in range(0, lat_size-1):
            
            div2[:, jj, :] = (1/np.sin(theta_rad[jj])) * (np.diff(u2[:, jj, :]) / diff_lamda_rad[::-1])  # 0<theta<pi

        #---------------------------------

        div3 = np.zeros((tim_size, lat_size-1, lon_size-1))
        
        for jj in range(0, lat_size-1):
            
            div3[:, jj, :] = (1/np.sin(theta_rad[jj])) * (div1[:, jj, :])   # 0<theta<pi

        #---------------------------------
            
        div4            = (1/r_sphere) * ( div2 + div3 ) 
        
        div4[:,  0,  :] = 0.5*( div4[:, -2,  :] + div4[:, 1, :] ) 
        
        div4[:,  :,  0] = 0.5*( div4[:,  :, -2] + div4[:, :, 1] )
        
        div4[:, -1,  :] = 0.5*( div4[:, -2,  :] + div4[:, 1, :] ) 
        
        div4[:,  :, -1] = 0.5*( div4[:,  :, -2] + div4[:, :, 1] )

        #--------------------------------- 
        
        data[exp]['div_layer_2_sum_glb'] = np.nansum(np.abs(div4[:,         :, :]), axis=(1,2))
        
        data[exp]['div_layer_2_plt_glb'] = (data[exp]['div_layer_2_sum_glb'])/1000	                # At the legend is added "(x 10^3)"
        
        #--------------------------------- 
    
        # Meridional mean 
        
        data[exp]['div_layer_2_mean_glb_inlat'] = np.nanmean(np.abs(div4[:, :, :]), axis=(0,1))*100	# At the legend is added "(x 10^-2)"
        
        # Zonal mean
        
        data[exp]['div_layer_2_mean_glb_inlon'] = np.nanmean(np.abs(div4[:, :, :]), axis=(0,2))*100	# At the legend is added "(x 10^-2)"
        
        #--------------------------------- 

        ##################################################################
        
    
    
    print()
    print('# Ploting -->  Sum of Divergence  =  Wave activity  (Layer 2)  -->  Simple line plot')
    print()
            
    fig, ax = plt.subplots(figsize=(10,5))

    plt.grid(True)
    
    plt.plot(tim_range, data['Dry_Baroclinic_Weak'  ]['div_layer_2_plt_glb'], '--', color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(tim_range, data['Dry_Baroclinic_Strong']['div_layer_2_plt_glb'], '--', color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Weak'  ]['div_layer_2_plt_glb'], '--', color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Strong']['div_layer_2_plt_glb'], '--', color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['div_layer_2_plt_glb'],       color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['div_layer_2_plt_glb'],       color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['div_layer_2_plt_glb'],       color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['div_layer_2_plt_glb'],       color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=12)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    ax.set_ylim(-0.2, 11)
    ax.yaxis.set_tick_params(labelsize=12)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.title(' Aeolus2.0 - Wave Activity in Layer 2 ', fontsize=15, pad=15)
    plt.ylabel('(x $10^{3}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Wave_Activity_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
        
        
        
    print()
    print('# Ploting -->  Meridional Mean of Divergence  (Layer 2)  -->  Simple line plot')
    print()
                
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)

    lon_plt = lon[:len(data['Dry_Baroclinic_Weak']['div_layer_2_mean_glb_inlat'])]

    plt.plot(lon_plt, data['Dry_Baroclinic_Weak'  ]['div_layer_2_mean_glb_inlat'], color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Baroclinic_Strong']['div_layer_2_mean_glb_inlat'], color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Weak'  ]['div_layer_2_mean_glb_inlat'], color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Strong']['div_layer_2_mean_glb_inlat'], color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Weak'   ]['div_layer_2_mean_glb_inlat'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Strong' ]['div_layer_2_mean_glb_inlat'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Weak'   ]['div_layer_2_mean_glb_inlat'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Strong' ]['div_layer_2_mean_glb_inlat'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    plt.xlabel('Longitude')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_xticks([-180, -135,  -90,  -45,    0,   45,   90,  135,  180])
    ax.set_xticklabels([u'180\N{DEGREE SIGN}', u'135\N{DEGREE SIGN}W', u'90\N{DEGREE SIGN}W', u'45\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'45\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E', u'135\N{DEGREE SIGN}E', u'180\N{DEGREE SIGN}'])
    ax.set_xlim(-180,180.05)
    
    xx = loc[0:2] # [-55, -35]
    y1, y2 = ax.get_ylim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.title(' Aeolus2.0 - Meridional mean of Divergence in Layer 2 ', fontsize=15, pad=15)
    plt.ylabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Divergence_Meridional_Mean_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    

    
    print()
    print('# Ploting -->  Zonal Mean of Divergence  (Layer 2)  -->  Simple line plot')
    print()
            
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)

    lat_plt = lat[:len(data['Dry_Baroclinic_Weak']['div_layer_2_mean_glb_inlon'])]

    plt.plot(data['Dry_Baroclinic_Weak'  ]['div_layer_2_mean_glb_inlon'], lat_plt, color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(data['Dry_Baroclinic_Strong']['div_layer_2_mean_glb_inlon'], lat_plt, color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(data['Dry_Barotropic_Weak'  ]['div_layer_2_mean_glb_inlon'], lat_plt, color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(data['Dry_Barotropic_Strong']['div_layer_2_mean_glb_inlon'], lat_plt, color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    
    plt.plot(data['MC_Baroclinic_Weak'   ]['div_layer_2_mean_glb_inlon'], lat_plt, color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(data['MC_Baroclinic_Strong' ]['div_layer_2_mean_glb_inlon'], lat_plt, color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(data['MC_Barotropic_Weak'   ]['div_layer_2_mean_glb_inlon'], lat_plt, color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(data['MC_Barotropic_Strong' ]['div_layer_2_mean_glb_inlon'], lat_plt, color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    plt.ylabel('Latitude')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
    ax.set_ylim(-90,90.05)
    
    y1 = [loc[2], loc[2]] # [10, 10]
    y2 = [loc[3], loc[3]] # [30, 30]
    xx = ax.get_xlim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')  
    
    plt.title(' Aeolus2.0 - Zonal mean of Divergence in Layer 2 ', fontsize=15, pad=15)
    plt.xlabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Divergence_Zonal_Mean_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)

    
    # End of plots_div_layer_2
    



plots_div_layer_1_2 = 1     # NEED to plots_div_layer_1=1 and plots_div_layer_2=1 FROM ABOVE


if plots_div_layer_1_2 == 1:
    
    
    print('')
    print('')
    print('> Divergence in Layer 1+2')
    print()

    
    print()
    print('# Ploting -->  Sum of Divergence  =  Wave activity  (Layer 1+2) --> Simple line plot')
    print()
            
    fig, ax = plt.subplots(figsize=(10,5))

    plt.grid(True)

    plt.plot(tim_range, (data['Dry_Baroclinic_Weak'  ]['div_layer_1_plt_glb'] + data['Dry_Baroclinic_Weak'  ]['div_layer_2_plt_glb']), '--', color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(tim_range, (data['Dry_Baroclinic_Strong']['div_layer_1_plt_glb'] + data['Dry_Baroclinic_Strong']['div_layer_2_plt_glb']), '--', color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(tim_range, (data['Dry_Barotropic_Weak'  ]['div_layer_1_plt_glb'] + data['Dry_Barotropic_Weak'  ]['div_layer_2_plt_glb']), '--', color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(tim_range, (data['Dry_Barotropic_Strong']['div_layer_1_plt_glb'] + data['Dry_Barotropic_Strong']['div_layer_2_plt_glb']), '--', color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(tim_range, (data['MC_Baroclinic_Weak'   ]['div_layer_1_plt_glb'] + data['MC_Baroclinic_Weak'   ]['div_layer_2_plt_glb']),       color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, (data['MC_Baroclinic_Strong' ]['div_layer_1_plt_glb'] + data['MC_Baroclinic_Strong' ]['div_layer_2_plt_glb']),       color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, (data['MC_Barotropic_Weak'   ]['div_layer_1_plt_glb'] + data['MC_Barotropic_Weak'   ]['div_layer_2_plt_glb']),       color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, (data['MC_Barotropic_Strong' ]['div_layer_1_plt_glb'] + data['MC_Barotropic_Strong' ]['div_layer_2_plt_glb']),       color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=12)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    ax.set_ylim(-0.3,15)
    ax.yaxis.set_tick_params(labelsize=12)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    
    plt.title(' Aeolus2.0 - Wave Activity in Layers 1+2 ', fontsize=15, pad=15)
    plt.ylabel('(x $10^{3}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Wave_Activity_Layers_1_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    
    
    # End of plots_div_layer_1_2




plots_hamiltonian = 1


if plots_hamiltonian == 1:
    
    
    print('')
    print('> Hamiltonian (Kinetic + Potential Energy) --> Layer 1 | Layer 2 | Layer 1 + Layer 2')
    print()
        

    for ie, exp in enumerate(cases):

        print('exp -->', exp)


        ##################################################################
        # Calculation of Hamiltonian Total Energy (Kinetic + Potential Energy)

        u1 = data[exp]['u1_t' ][:,::-1,:] # u1_t[:,:,:]     # u1ph
        v1 = data[exp]['v1_t' ][:,::-1,:] # v1_t[:,:,:]     # u1th

        u2 = data[exp]['u2_t' ][:,::-1,:] # u2_t[:,:,:]     # u2ph
        v2 = data[exp]['v2_t' ][:,::-1,:] # v2_t[:,:,:]     # u2th

        h1 = data[exp]['h1_t' ][:,::-1,:] # h1_t[:,:,:]     
        h2 = data[exp]['h2_t' ][:,::-1,:] # h2_t[:,:,:]     

        b1 = data[exp]['b1_t' ][:,::-1,:] # b1_t[:,:,:]     
        b2 = data[exp]['b2_t' ][:,::-1,:] # b2_t[:,:,:]     
        
        # The -1 above is to reverse the matrix to conform with matlab version
        
        # This was checked with Masoud by the paper 
        #--------------------------------- 
            
        h1_tild                 = 0.5*(h1+H1) + (h2+H2)
        
        energy_layer_1          = (h1+H1)*( 0.5*(u1**2 + v1**2) + h1_tild*(b1+B1) )

        energy_layer_1_sum      = np.nansum(energy_layer_1, axis=0)

        #---------------------------------
        
        h2_tild                 = 0.5*(h2+H2)
            
        energy_layer_2          = (h2+H2)*( 0.5*(u2**2 + v2**2) + h2_tild*(b2+B2) ) 

        energy_layer_2_sum      = np.nansum(energy_layer_2, axis=0)
                
        #--------------------------------- 
        
        energy_layer_all        = energy_layer_1 + energy_layer_2
        
        energy_layer_all_timsum = np.nansum(energy_layer_all, axis=0) 
        
        #---------------------------------
                
        if ie == 0:
            
            data[exp]['energy_layer_all_fldsum_glb'     ] = np.nansum(energy_layer_all[:, :, :], axis=(1,2))
            data[exp]['energy_layer_all_fldsum_glb_perc'] =   ( data[exp]['energy_layer_all_fldsum_glb'] - data[exp]['energy_layer_all_fldsum_glb'][0] ) / data[exp]['energy_layer_all_fldsum_glb'][0]
            
            data[exp]['energy_layer_1_fldsum_glb'       ] = np.nansum(energy_layer_1[:, :, :],   axis=(1,2))
            data[exp]['energy_layer_1_fldsum_glb_perc'  ] =   ( data[exp]['energy_layer_1_fldsum_glb']   - data[exp]['energy_layer_1_fldsum_glb'][0] )   / data[exp]['energy_layer_1_fldsum_glb'][0]
            
            data[exp]['energy_layer_2_fldsum_glb'       ] = np.nansum(energy_layer_2[:, :, :],   axis=(1,2))
            data[exp]['energy_layer_2_fldsum_glb_perc'  ] =   ( data[exp]['energy_layer_2_fldsum_glb']   - data[exp]['energy_layer_2_fldsum_glb'][0] )   / data[exp]['energy_layer_2_fldsum_glb'][0]
            
            CTE_energy_layer_all_fldsum_glb      = data[exp]['energy_layer_all_fldsum_glb'     ][0]
            CTE_energy_layer_all_fldsum_glb_perc = data[exp]['energy_layer_all_fldsum_glb_perc'][0]
            
            CTE_energy_layer_1_fldsum_glb        = data[exp]['energy_layer_1_fldsum_glb'     ][0]
            CTE_energy_layer_1_fldsum_glb_perc   = data[exp]['energy_layer_1_fldsum_glb_perc'][0]
            
            CTE_energy_layer_2_fldsum_glb        = data[exp]['energy_layer_2_fldsum_glb'     ][0]
            CTE_energy_layer_2_fldsum_glb_perc   = data[exp]['energy_layer_2_fldsum_glb_perc'][0]


        data[exp]['energy_layer_all_fldsum_glb'     ] = np.nansum(energy_layer_all[:, :, :], axis=(1,2))
        data[exp]['energy_layer_all_fldsum_glb_perc'] = ( ( data[exp]['energy_layer_all_fldsum_glb'] - CTE_energy_layer_all_fldsum_glb ) / CTE_energy_layer_all_fldsum_glb )*1000 # + ( CTE_energy_layer_all_fldsum_glb_perc - ( data[exp]['energy_layer_all_fldsum_glb'][0] /CTE_energy_layer_all_fldsum_glb) ) ) *100

        data[exp]['energy_layer_1_fldsum_glb'       ] = np.nansum(energy_layer_1[:, :, :],   axis=(1,2))
        data[exp]['energy_layer_1_fldsum_glb_perc'  ] = ( ( data[exp]['energy_layer_1_fldsum_glb'  ] - CTE_energy_layer_1_fldsum_glb   ) / CTE_energy_layer_all_fldsum_glb )*1000 # + ( CTE_energy_layer_all_fldsum_glb_perc - ( data[exp]['energy_layer_1_fldsum_glb'  ][0] /CTE_energy_layer_all_fldsum_glb) ) ) *100

        data[exp]['energy_layer_2_fldsum_glb'       ] = np.nansum(energy_layer_2[:, :, :],   axis=(1,2))
        data[exp]['energy_layer_2_fldsum_glb_perc'  ] = ( ( data[exp]['energy_layer_2_fldsum_glb'  ] - CTE_energy_layer_2_fldsum_glb   ) / CTE_energy_layer_all_fldsum_glb )*1000 # + ( CTE_energy_layer_all_fldsum_glb_perc - ( data[exp]['energy_layer_2_fldsum_glb'  ][0] /CTE_energy_layer_all_fldsum_glb) ) ) *100

        #---------------------------------
                
        ##################################################################
                
        
        
    print()
    print('# Ploting -->  Sum of Energy  =  Hamiltonian Total Energy (Layer 1 + Layer 2) --> Simple line plot --> time x Energy values')
    print()

    fig, ax = plt.subplots(figsize=(10,5))

    plt.grid(True)

    plt.plot(tim_range, data['Dry_Baroclinic_Weak'  ]['energy_layer_all_fldsum_glb_perc'], '--', color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(tim_range, data['Dry_Baroclinic_Strong']['energy_layer_all_fldsum_glb_perc'], '--', color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Weak'  ]['energy_layer_all_fldsum_glb_perc'], '--', color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Strong']['energy_layer_all_fldsum_glb_perc'], '--', color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['energy_layer_all_fldsum_glb_perc'],       color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['energy_layer_all_fldsum_glb_perc'],       color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['energy_layer_all_fldsum_glb_perc'],       color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['energy_layer_all_fldsum_glb_perc'],       color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=12)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    ax.set_ylim(-1, 5)
    ax.yaxis.set_tick_params(labelsize=12)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.ylabel('Hamiltonian Anomaly [$(H-H_{t=0})/H_{t=0}$] (‰)', fontsize=14)
    
    plt.title(' Aeolus2.0 - Hamiltonian (Kinetic + Potential Energy) in Layers 1+2 ', fontsize=15, pad=15)  
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Hamiltonian_Energy_Layers_1_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    
    
    
    print()
    print('# Ploting -->  Sum of Energy  =  Hamiltonian Total Energy (Layer 1) --> Simple line plot --> time x Energy values')
    print()

    fig, ax = plt.subplots(figsize=(10,5))

    plt.grid(True)

    plt.plot(tim_range, data['Dry_Baroclinic_Weak'  ]['energy_layer_1_fldsum_glb_perc'], '--', color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(tim_range, data['Dry_Baroclinic_Strong']['energy_layer_1_fldsum_glb_perc'], '--', color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Weak'  ]['energy_layer_1_fldsum_glb_perc'], '--', color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Strong']['energy_layer_1_fldsum_glb_perc'], '--', color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['energy_layer_1_fldsum_glb_perc'],       color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['energy_layer_1_fldsum_glb_perc'],       color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['energy_layer_1_fldsum_glb_perc'],       color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['energy_layer_1_fldsum_glb_perc'],       color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=12)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    ax.set_ylim(-1, 5)
    ax.yaxis.set_tick_params(labelsize=12)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.ylabel('Hamiltonian Anomaly [$(H-H_{t=0})/H_{t=0}$] (‰)', fontsize=14)
    
    plt.title(' Aeolus2.0 - Hamiltonian (Kinetic + Potential Energy) in Layer 1 ', fontsize=15, pad=15) 
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Hamiltonian_Energy_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)


    
    print()
    print('# Ploting -->  Sum of Energy  =  Hamiltonian Total Energy (Layer 2) --> Simple line plot --> time x Energy values')
    print()

    fig, ax = plt.subplots(figsize=(10,5))

    plt.grid(True)

    plt.plot(tim_range, data['Dry_Baroclinic_Weak'  ]['energy_layer_2_fldsum_glb_perc'], '--', color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(tim_range, data['Dry_Baroclinic_Strong']['energy_layer_2_fldsum_glb_perc'], '--', color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Weak'  ]['energy_layer_2_fldsum_glb_perc'], '--', color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(tim_range, data['Dry_Barotropic_Strong']['energy_layer_2_fldsum_glb_perc'], '--', color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['energy_layer_2_fldsum_glb_perc'],       color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['energy_layer_2_fldsum_glb_perc'],       color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['energy_layer_2_fldsum_glb_perc'],       color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['energy_layer_2_fldsum_glb_perc'],       color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=12)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    ax.set_ylim(-1, 5)
    ax.yaxis.set_tick_params(labelsize=12)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.ylabel('Hamiltonian Anomaly [$(H-H_{t=0})/H_{t=0}$] (‰)', fontsize=14)
    
    plt.title(' Aeolus2.0 - Hamiltonian (Kinetic + Potential Energy) in Layer 2 ', fontsize=15, pad=15) 
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_Hamiltonian_Energy_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    
    
    # End of plots_hamiltonian




plots_wind_layer_1 = 1


if plots_wind_layer_1 == 1:
    
    
    print('')
    print('> Wind U and V in Layer 1')
    print()
        

    for ie, exp in enumerate(cases):

        print('exp -->', exp)

        ##################################################################
        # Wind

        u1 = data[exp]['u1_t' ] # u1_t[:,:,:]   # u1ph
        v1 = data[exp]['v1_t' ] # v1_t[:,:,:]   # u1th

        #--------------------------------- 
        # Meridional
        
        data[exp]['u1_mean_glb_inlat'] = np.nanmean(u1[:, :, :], axis=(0,1))*100	# At the legend is added "(x 10^-2)"
        
        data[exp]['v1_mean_glb_inlat'] = np.nanmean(v1[:, :, :], axis=(0,1))*100	# At the legend is added "(x 10^-2)"
        
        #--------------------------------- 
        # Zonal
        
        data[exp]['u1_mean_glb_inlon'] = np.nanmean(u1[:, :, :], axis=(0,2))*100	# At the legend is added "(x 10^-2)"
        
        data[exp]['v1_mean_glb_inlon'] = np.nanmean(v1[:, :, :], axis=(0,2))*100	# At the legend is added "(x 10^-2)"
        
        ##################################################################
        
    
    
    print()
    print('# Ploting -->  Meridional mean of U (Layer 1) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(8,7))

    plt.grid(True)
    
    lon_plt = lon[:len(data['Dry_Baroclinic_Weak']['u1_mean_glb_inlat'])]

    plt.plot(lon_plt, data['Dry_Baroclinic_Weak'  ]['u1_mean_glb_inlat'], color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Baroclinic_Strong']['u1_mean_glb_inlat'], color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Weak'  ]['u1_mean_glb_inlat'], color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Strong']['u1_mean_glb_inlat'], color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Weak'   ]['u1_mean_glb_inlat'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Strong' ]['u1_mean_glb_inlat'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Weak'   ]['u1_mean_glb_inlat'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Strong' ]['u1_mean_glb_inlat'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11, loc='lower right')

    ax.set_xticks([-180, -135,  -90,  -45,    0,   45,   90,  135,  180])
    ax.set_xticklabels([u'180\N{DEGREE SIGN}', u'135\N{DEGREE SIGN}W', u'90\N{DEGREE SIGN}W', u'45\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'45\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E', u'135\N{DEGREE SIGN}E', u'180\N{DEGREE SIGN}'])
    ax.set_xlim(-180,180.05)
    ax.set_ylim(-0.17,0.15)

    xx = loc[0:2] # [-55, -35]
    y1, y2 = ax.get_ylim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.xlabel('Longitude')

    plt.title(' Aeolus2.0 - Meridional mean of U wind in Layer 1', fontsize=15, pad=15)
    plt.ylabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_U_Meridional_Mean_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)
    
    
    
    print()
    print('# Ploting -->  Meridional mean of V (Layer 1) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)
    
    lon_plt = lon[:len(data['Dry_Baroclinic_Weak']['v1_mean_glb_inlat'])]

    plt.plot(lon_plt, data['Dry_Baroclinic_Weak'  ]['v1_mean_glb_inlat'], color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Baroclinic_Strong']['v1_mean_glb_inlat'], color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Weak'  ]['v1_mean_glb_inlat'], color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Strong']['v1_mean_glb_inlat'], color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Weak'   ]['v1_mean_glb_inlat'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Strong' ]['v1_mean_glb_inlat'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Weak'   ]['v1_mean_glb_inlat'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Strong' ]['v1_mean_glb_inlat'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_xticks([-180, -135,  -90,  -45,    0,   45,   90,  135,  180])
    ax.set_xticklabels([u'180\N{DEGREE SIGN}', u'135\N{DEGREE SIGN}W', u'90\N{DEGREE SIGN}W', u'45\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'45\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E', u'135\N{DEGREE SIGN}E', u'180\N{DEGREE SIGN}'])
    ax.set_xlim(-180,180.05)
    ax.set_ylim(-0.12,0.25)

    xx = loc[0:2] # [-55, -35]
    y1, y2 = ax.get_ylim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.xlabel('Longitude')

    plt.title(' Aeolus2.0 - Meridional mean of V wind in Layer 1', fontsize=15, pad=15)
    plt.ylabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_V_Meridional_Mean_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)



    
    print()
    print('# Ploting -->  Zonal mean of U (Layer 1) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)
    
    lat_plt = lat[:len(data['Dry_Baroclinic_Weak']['u1_mean_glb_inlon'])]
    
    plt.plot(data['Dry_Baroclinic_Weak'  ]['u1_mean_glb_inlon'], lat_plt, color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(data['Dry_Baroclinic_Strong']['u1_mean_glb_inlon'], lat_plt, color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(data['Dry_Barotropic_Weak'  ]['u1_mean_glb_inlon'], lat_plt, color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(data['Dry_Barotropic_Strong']['u1_mean_glb_inlon'], lat_plt, color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    
    plt.plot(data['MC_Baroclinic_Weak'   ]['u1_mean_glb_inlon'], lat_plt, color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(data['MC_Baroclinic_Strong' ]['u1_mean_glb_inlon'], lat_plt, color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(data['MC_Barotropic_Weak'   ]['u1_mean_glb_inlon'], lat_plt, color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(data['MC_Barotropic_Strong' ]['u1_mean_glb_inlon'], lat_plt, color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
    ax.set_ylim(-90,90.05)
    ax.set_xlim(-0.5,0.5)

    y1 = [loc[2], loc[2]] # [10, 10]
    y2 = [loc[3], loc[3]] # [30, 30]
    xx = ax.get_xlim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.ylabel('Latitude')

    plt.title(' Aeolus2.0 - Zonal mean of U wind in Layer 1', fontsize=15, pad=15)
    plt.xlabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_U_Zonal_Mean_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)


    
    print()
    print('# Ploting -->  Zonal mean of V (Layer 1) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(8,7))

    plt.grid(True)
    
    lat_plt = lat[:len(data['Dry_Baroclinic_Weak']['v1_mean_glb_inlon'])]
    
    plt.plot(data['Dry_Baroclinic_Weak'  ]['v1_mean_glb_inlon'], lat_plt, color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(data['Dry_Baroclinic_Strong']['v1_mean_glb_inlon'], lat_plt, color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(data['Dry_Barotropic_Weak'  ]['v1_mean_glb_inlon'], lat_plt, color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(data['Dry_Barotropic_Strong']['v1_mean_glb_inlon'], lat_plt, color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    
    plt.plot(data['MC_Baroclinic_Weak'   ]['v1_mean_glb_inlon'], lat_plt, color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(data['MC_Baroclinic_Strong' ]['v1_mean_glb_inlon'], lat_plt, color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(data['MC_Barotropic_Weak'   ]['v1_mean_glb_inlon'], lat_plt, color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(data['MC_Barotropic_Strong' ]['v1_mean_glb_inlon'], lat_plt, color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11, loc='lower left')

    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
    ax.set_ylim(-90,90.05)
    ax.set_xlim(-0.03,0.03)

    y1 = [loc[2], loc[2]] # [10, 10]
    y2 = [loc[3], loc[3]] # [30, 30]
    xx = ax.get_xlim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.ylabel('Latitude')

    plt.title(' Aeolus2.0 - Zonal mean of V wind in Layer 1', fontsize=15, pad=15)
    plt.xlabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_V_Zonal_Mean_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)


    # End of plots_wind_layer_1




plots_wind_layer_2 = 1


if plots_wind_layer_2 == 1:
    
    
    print('')
    print('> Wind U and V in Layer 2')
    print()
        

    for ie, exp in enumerate(cases):

        print('exp -->', exp)

        ##################################################################
        # Wind

        u2 = data[exp]['u2_t' ] # u2_t[:,:,:]   # u2ph
        v2 = data[exp]['v2_t' ] # v2_t[:,:,:]   # u2th

        #--------------------------------- 
        # Meridional
        
        data[exp]['u2_mean_glb_inlat'] = np.nanmean(u2[:, :, :], axis=(0,1))*100	# At the legend is added "(x 10^-2)"
        
        data[exp]['v2_mean_glb_inlat'] = np.nanmean(v2[:, :, :], axis=(0,1))*100	# At the legend is added "(x 10^-2)"
        
        #--------------------------------- 
        # Zonal
        
        data[exp]['u2_mean_glb_inlon'] = np.nanmean(u2[:, :, :], axis=(0,2))*100	# At the legend is added "(x 10^-2)"
        
        data[exp]['v2_mean_glb_inlon'] = np.nanmean(v2[:, :, :], axis=(0,2))*100	# At the legend is added "(x 10^-2)"
        
        ##################################################################
        
    
    
    print()
    print('# Ploting -->  Meridional mean of U (Layer 2) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(8,7))

    plt.grid(True)
    
    lon_plt = lon[:len(data['Dry_Baroclinic_Weak']['u2_mean_glb_inlat'])]

    plt.plot(lon_plt, data['Dry_Baroclinic_Weak'  ]['u2_mean_glb_inlat'], color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Baroclinic_Strong']['u2_mean_glb_inlat'], color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Weak'  ]['u2_mean_glb_inlat'], color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Strong']['u2_mean_glb_inlat'], color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Weak'   ]['u2_mean_glb_inlat'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Strong' ]['u2_mean_glb_inlat'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Weak'   ]['u2_mean_glb_inlat'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Strong' ]['u2_mean_glb_inlat'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11, loc='lower right')

    ax.set_xticks([-180, -135,  -90,  -45,    0,   45,   90,  135,  180])
    ax.set_xticklabels([u'180\N{DEGREE SIGN}', u'135\N{DEGREE SIGN}W', u'90\N{DEGREE SIGN}W', u'45\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'45\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E', u'135\N{DEGREE SIGN}E', u'180\N{DEGREE SIGN}'])
    ax.set_xlim(-180,180.05)
    ax.set_ylim(-0.17,0.15)

    xx = loc[0:2] # [-55, -35]
    y1, y2 = ax.get_ylim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.xlabel('Longitude')

    plt.title(' Aeolus2.0 - Meridional mean of U wind in Layer 2', fontsize=15, pad=15)
    plt.ylabel('(x $10^{-2}$)')
        
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_U_Meridional_Mean_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)
    
    
    
    print()
    print('# Ploting -->  Meridional mean of V (Layer 2) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)
    
    lon_plt = lon[:len(data['Dry_Baroclinic_Weak']['v2_mean_glb_inlat'])]

    plt.plot(lon_plt, data['Dry_Baroclinic_Weak'  ]['v2_mean_glb_inlat'], color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Baroclinic_Strong']['v2_mean_glb_inlat'], color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Weak'  ]['v2_mean_glb_inlat'], color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(lon_plt, data['Dry_Barotropic_Strong']['v2_mean_glb_inlat'], color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Weak'   ]['v2_mean_glb_inlat'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(lon_plt, data['MC_Baroclinic_Strong' ]['v2_mean_glb_inlat'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Weak'   ]['v2_mean_glb_inlat'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(lon_plt, data['MC_Barotropic_Strong' ]['v2_mean_glb_inlat'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_xticks([-180, -135,  -90,  -45,    0,   45,   90,  135,  180])
    ax.set_xticklabels([u'180\N{DEGREE SIGN}', u'135\N{DEGREE SIGN}W', u'90\N{DEGREE SIGN}W', u'45\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'45\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E', u'135\N{DEGREE SIGN}E', u'180\N{DEGREE SIGN}'])
    ax.set_xlim(-180,180.05)
    ax.set_ylim(-0.35,0.2)

    xx = loc[0:2] # [-55, -35]
    y1, y2 = ax.get_ylim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.xlabel('Longitude')

    plt.title(' Aeolus2.0 - Meridional mean of V wind in Layer 2', fontsize=15, pad=15)
    plt.ylabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_V_Meridional_Mean_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)



    
    print()
    print('# Ploting -->  Zonal mean of U (Layer 2) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(10,9))

    plt.grid(True)
    
    lat_plt = lat[:len(data['Dry_Baroclinic_Weak']['u2_mean_glb_inlon'])]
    
    plt.plot(data['Dry_Baroclinic_Weak'  ]['u2_mean_glb_inlon'], lat_plt, color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(data['Dry_Baroclinic_Strong']['u2_mean_glb_inlon'], lat_plt, color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(data['Dry_Barotropic_Weak'  ]['u2_mean_glb_inlon'], lat_plt, color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(data['Dry_Barotropic_Strong']['u2_mean_glb_inlon'], lat_plt, color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    
    plt.plot(data['MC_Baroclinic_Weak'   ]['u2_mean_glb_inlon'], lat_plt, color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(data['MC_Baroclinic_Strong' ]['u2_mean_glb_inlon'], lat_plt, color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(data['MC_Barotropic_Weak'   ]['u2_mean_glb_inlon'], lat_plt, color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(data['MC_Barotropic_Strong' ]['u2_mean_glb_inlon'], lat_plt, color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)

    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
    ax.set_ylim(-90,90.05)
    ax.set_xlim(-0.5,0.5)

    y1 = [loc[2], loc[2]] # [10, 10]
    y2 = [loc[3], loc[3]] # [30, 30]
    xx = ax.get_xlim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.ylabel('Latitude')

    plt.title(' Aeolus2.0 - Zonal mean of U wind in Layer 2', fontsize=15, pad=15)
    plt.xlabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_U_Zonal_Mean_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)


    
    print()
    print('# Ploting -->  Zonal mean of V (Layer 2) --> Simple line plot --> lon x wind values')
    print()
    
    fig, ax = plt.subplots(figsize=(8,7))

    plt.grid(True)
    
    lat_plt = lat[:len(data['Dry_Baroclinic_Weak']['v2_mean_glb_inlon'])]
    
    plt.plot(data['Dry_Baroclinic_Weak'  ]['v2_mean_glb_inlon'], lat_plt, color='yellow',        linewidth=3, label='Dry Baroclinic Weak   Bump')
    plt.plot(data['Dry_Baroclinic_Strong']['v2_mean_glb_inlon'], lat_plt, color='orange',        linewidth=3, label='Dry Baroclinic Strong Bump')
    plt.plot(data['Dry_Barotropic_Weak'  ]['v2_mean_glb_inlon'], lat_plt, color='palevioletred', linewidth=3, label='Dry Barotropic Weak   Bump')
    plt.plot(data['Dry_Barotropic_Strong']['v2_mean_glb_inlon'], lat_plt, color='red',           linewidth=3, label='Dry Barotropic Strong Bump')
    
    plt.plot(data['MC_Baroclinic_Weak'   ]['v2_mean_glb_inlon'], lat_plt, color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(data['MC_Baroclinic_Strong' ]['v2_mean_glb_inlon'], lat_plt, color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(data['MC_Barotropic_Weak'   ]['v2_mean_glb_inlon'], lat_plt, color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(data['MC_Barotropic_Strong' ]['v2_mean_glb_inlon'], lat_plt, color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11, loc='lower left')

    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
    ax.set_ylim(-90,90.05)
    ax.set_xlim(-0.1,0.03)

    y1 = [loc[2], loc[2]] # [10, 10]
    y2 = [loc[3], loc[3]] # [30, 30]
    xx = ax.get_xlim()
    ax.fill_between(xx, y1, y2, alpha=0.95, color='lightgray')
    
    plt.ylabel('Latitude')

    plt.title(' Aeolus2.0 - Zonal mean of V wind in Layer 2 ', fontsize=15, pad=15)
    plt.xlabel('(x $10^{-2}$)')
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Plot_V_Zonal_Mean_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

    if exists(figname): print('done -->', figname)


    # End of plots_wind_layer_2




plots_clwc_layer_1 = 1


if plots_clwc_layer_1 == 1:


    print('')
    print('> CLWC in Layer 1')
    print()
        

    for ie, exp in enumerate(cases):

        print('exp -->', exp)

        ##################################################################
        
        clwc = data[exp]['cc1_t']  # Condensed liquid water content (CLWC)
        
        data[exp]['clwc_fldsum'] = np.nansum(clwc, axis=(1,2))
        
        ##################################################################
        
    
    print()
    print('# Ploting -->  Spatial sum of CLWC (Layer 1) --> Simple line plot --> time x values')
    print()
    
    fig, ax = plt.subplots(figsize=(10,6))

    plt.grid(True)
        
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['clwc_fldsum'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['clwc_fldsum'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['clwc_fldsum'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['clwc_fldsum'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=8)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    # ax.set_ylim(0.0,1.05)
    ylim = ax.get_ylim()
    
    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.title(' Aeolus2.0 - Condensed Liquid Water Content (CLWC) in Layer 1 ', fontsize=15, pad=15)
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Water_CLWC_Layer_1.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    
    
    # End of plots_clwc_layer_1


    

plots_w2_layer_2 = 1


if plots_w2_layer_2 == 1:
    
    
    print('')
    print('> W2 in Layer 2')
    print()
        

    for ie, exp in enumerate(cases):

        print('exp -->', exp)

        ##################################################################
        
        w2 = data[exp]['w2_t']  # Bulk of precipitable water
        
        data[exp]['w2_fldsum'] = np.nansum(w2, axis=(1,2))
        
        ##################################################################
        
    
    print()
    print('# Ploting -->  Spatial sum of W2 (Layer 2) --> Simple line plot --> time x values')
    print()
    
    fig, ax = plt.subplots(figsize=(10,6))

    plt.grid(True)
        
    plt.plot(tim_range, data['MC_Baroclinic_Weak'   ]['w2_fldsum'], color='deepskyblue',   linewidth=1, label='MC Baroclinic Weak   Bump')
    plt.plot(tim_range, data['MC_Baroclinic_Strong' ]['w2_fldsum'], color='blue',          linewidth=1, label='MC Baroclinic Strong Bump')
    plt.plot(tim_range, data['MC_Barotropic_Weak'   ]['w2_fldsum'], color='gray',          linewidth=1, label='MC Barotropic Weak   Bump')
    plt.plot(tim_range, data['MC_Barotropic_Strong' ]['w2_fldsum'], color='black',         linewidth=1, label='MC Barotropic Strong Bump')

    ax.legend(fancybox=True, framealpha=0.2, fontsize=11)
    ax.set_xticks( tim_range[19::20])
    ax.set_xticklabels(times[19::20])
    ax.xaxis.set_tick_params(labelsize=8)
    plt.xticks(rotation=45, ha='right')
    ax.set_xlim(tim_range[0], tim_range[-1]+10)
    # ax.set_ylim(0.0,1.05)
    ylim = ax.get_ylim()

    plt.plot([tim_range[21 ], tim_range[21 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[26 ], tim_range[26 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[35 ], tim_range[35 ]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[111], tim_range[111]], ylim, '--', color='black', alpha=0.3, linewidth=2)
    plt.plot([tim_range[199], tim_range[199]], ylim, '--', color='black', alpha=0.3, linewidth=2)

    plt.title(' Aeolus2.0 - Bulk of Precipitable Water (W2) in Layer 2 ', fontsize=15, pad=15)
    
    # Saving #########################################################

    figname = '{}/Aeolus2.0_Water_W2_Layer_2.png'.format(dir_figs)

    fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

    if exists(figname): print('done -->', figname)
    
            
    # End of plots_w2_layer_2










print('')
print('')
print('')
print('')
print('')
print('The End')
print('')
print('')


