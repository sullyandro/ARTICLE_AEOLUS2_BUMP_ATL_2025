################################################################################
#
# Author:       Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Colaborators: Masoud Rostami
# Date:         01.08.2025
# Type:         Python3
#
# Description:
# Script to plot Hovmoller graphics for Aeolus2 analysis of Bump evolution.
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
print('Aeolus2.0 Analysis Preparation - Plot Hovmoller')


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


for exp in cases:

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

        if var == 'u1ph': u1_t  = dat.variables['u1ph'][:] # Zonal (azimuthal) velocity layer 1           ; units: [(g*H)^0.5], [beta*L_d^2] at the equator
        if var == 'u1th': v1_t  = dat.variables['u1th'][:] # Meridional velocity layer 1                  ; units: [(g*H)^0.5], [beta*L_d^2] at the equator

        if var == 'u2ph': u2_t  = dat.variables['u2ph'][:] # Zonal (azimuthal) velocity layer 2           ; units: [(g*H)^0.5], [beta*L_d^2] at the equator
        if var == 'u2th': v2_t  = dat.variables['u2th'][:] # Meridional velocity layer 2                  ; units: [(g*H)^0.5], [beta*L_d^2] at the equator

        if var == 'h1':   h1_t  = dat.variables['h1'  ][:] # Pseudo-height layer 1                        ; units: H
        if var == 'h2':   h2_t  = dat.variables['h2'  ][:] # Pseudo-height layer 2                        ; units: H

        if var == 'b1':   b1_t  = dat.variables['b1'  ][:] # Buoyancy layer 1                             ; units: g*theta/theta_s
        if var == 'b2':   b2_t  = dat.variables['b2'  ][:] # Buoyancy layer 2                             ; units: g*theta/theta_s

        if var == 'CC1':  cc1_t = dat.variables['CC1' ][:] # CLWC (Cloud Liquid Water Content) at layer 1 ; units: Q/T
        if var == 'w2':   w2_t  = dat.variables['w2'  ][:] # Bulk of Precipitable Water at layer 2        ; units: (L_v.g/(C_p.theta_s))Kg/Kg

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
    print('# Plots Hovmoller production')
    print()




    if remap == 0: dir_figs = '../Figures/Plot_Hovmoller/Aeolus2.0_{}'.format(exp)
    if remap == 1: dir_figs = '../Figures_remapbil_1dg/Plot_Hovmoller/Aeolus2.0_{}'.format(exp)

    if not exists(dir_figs): os.makedirs(dir_figs)

    
    
    
        
    plots_div_layer_1 = 1


    if plots_div_layer_1 == 1:
        
        print('')
        print('> Divergence in Layer 1')
        print()
        

        ##################################################################
        # Calculation of Divergence

        u1 = u1_t[:,:,:]    # u1ph
        v1 = v1_t[:,:,:]    # u1th

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
        
        # Meridional mean 
        
        div_layer_1_mean_glb_inlat      = np.nanmean(np.abs(div4[:, :, :]), axis=(0,1))*100 # At the legend is added "(x 10^-2)"
        div_layer_1_mean_glb_inlat_hov  = np.nanmean(np.abs(div4[:, :, :]), axis=1    )*100 # At the legend is added "(x 10^-2)"
        
        # Zonal mean
        
        div_layer_1_mean_glb_inlon      = np.nanmean(np.abs(div4[:, :, :]), axis=(0,2))*100 # At the legend is added "(x 10^-2)"
        div_layer_1_mean_glb_inlon_hov  = np.nanmean(np.abs(div4[:, :, :]), axis=2    )*100 # At the legend is added "(x 10^-2)"
        
        #--------------------------------- 
        
        ##################################################################


        
        print()
        print('# Ploting -->  Meridional Mean of Divergence  (Layer 1)  -->  Hovmöller plot  -->  time x lon - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(0.2, 2, 0.2)
        clevs = np.concatenate([[0], clevs])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lon[:len(div_layer_1_mean_glb_inlat_hov[0,:])], np.swapaxes(div_layer_1_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='max')

        cs = ax2.contour (dtime_hov, lon[:len(div_layer_1_mean_glb_inlat_hov[0,:])], np.swapaxes(div_layer_1_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
        ax2.set_yticklabels([u'', u'120\N{DEGREE SIGN}W', u'60\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'60\N{DEGREE SIGN}E', u'120\N{DEGREE SIGN}E', u''][::-1])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[0], -loc[0]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[1], -loc[1]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)

        plt.title(' Aeolus2.0 - Case {} \n Meridional mean of Divergence in Layer 1'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Divergence_Meridional_Mean_Layer_1.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

        if exists(figname): print('done -->', figname)
                
        

        print()
        print('# Ploting -->  Zonal Mean of Divergence  (Layer 1)  -->  Hovmöller plot  -->  time x lat - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(0.2, 2, 0.2)
        clevs = np.concatenate([[0], clevs])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lat[:len(div_layer_1_mean_glb_inlon_hov[0,:])], np.swapaxes(div_layer_1_mean_glb_inlon_hov, 0,1), clevs, cmap=cmap_new, extend='max')

        cs = ax2.contour (dtime_hov, lat[:len(div_layer_1_mean_glb_inlon_hov[0,:])], np.swapaxes(div_layer_1_mean_glb_inlon_hov, 0,1), clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax2.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[2], loc[2]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[3], loc[3]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.title(' Aeolus2.0 - Case {} \n Zonal mean of Divergence in Layer 1'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Divergence_Zonal_Mean_Layer_1.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.close()

        if exists(figname): print('done -->', figname)
        
        
    # End of plots_div_layer_1
    



    plots_div_layer_2 = 1


    if plots_div_layer_2 == 1:
        
        print('')
        print('> Divergence in Layer 2')
        print()
        

        ##################################################################
        # Calculation of Divergence

        u2 = u2_t[:,:,:]    # u2ph
        v2 = v2_t[:,:,:]    # u2th

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
        
        # Meridional mean 
        
        div_layer_2_mean_glb_inlat      = np.nanmean(np.abs(div4[:, :, :]), axis=(0,1))*100 # At the legend is added "(x 10^-2)"
        div_layer_2_mean_glb_inlat_hov  = np.nanmean(np.abs(div4[:, :, :]), axis=1    )*100 # At the legend is added "(x 10^-2)"
        
        # Zonal mean
        
        div_layer_2_mean_glb_inlon      = np.nanmean(np.abs(div4[:, :, :]), axis=(0,2))*100 # At the legend is added "(x 10^-2)"
        div_layer_2_mean_glb_inlon_hov  = np.nanmean(np.abs(div4[:, :, :]), axis=2    )*100 # At the legend is added "(x 10^-2)"
        
        #--------------------------------- 
        
        ##################################################################


        
        print()
        print('# Ploting -->  Meridional Mean of Divergence  (Layer 2)  -->  Hovmöller plot  -->  time x lon - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(0.2, 2, 0.2)
        clevs = np.concatenate([[0], clevs])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lon[:len(div_layer_2_mean_glb_inlat_hov[0,:])], np.swapaxes(div_layer_2_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='max')

        cs = ax2.contour (dtime_hov, lon[:len(div_layer_2_mean_glb_inlat_hov[0,:])], np.swapaxes(div_layer_2_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
        ax2.set_yticklabels([u'', u'120\N{DEGREE SIGN}W', u'60\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'60\N{DEGREE SIGN}E', u'120\N{DEGREE SIGN}E', u''][::-1])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[0], -loc[0]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[1], -loc[1]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)

        plt.title(' Aeolus2.0 - Case {} \n Meridional mean of Divergence in Layer 2'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Divergence_Meridional_Mean_Layer_2.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

        if exists(figname): print('done -->', figname)
                
        

        print()
        print('# Ploting -->  Zonal Mean of Divergence  (Layer 2)  -->  Hovmöller plot  -->  time x lat - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(0.2, 2, 0.2)
        clevs = np.concatenate([[0], clevs])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lat[:len(div_layer_2_mean_glb_inlon_hov[0,:])], np.swapaxes(div_layer_2_mean_glb_inlon_hov, 0,1), clevs, cmap=cmap_new, extend='max')

        cs = ax2.contour (dtime_hov, lat[:len(div_layer_2_mean_glb_inlon_hov[0,:])], np.swapaxes(div_layer_2_mean_glb_inlon_hov, 0,1), clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax2.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[2], loc[2]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[3], loc[3]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.title(' Aeolus2.0 - Case {} \n Zonal mean of Divergence in Layer 2'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Divergence_Zonal_Mean_Layer_2.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.close()

        if exists(figname): print('done -->', figname)
        
        
    # End of plots_div_layer_2
    
    


    plots_hamiltonian_layer_1 = 1


    if plots_hamiltonian_layer_1 == 1:
        
        print('')
        print('> Hamiltonian (Kinetic + Potential Energy) in Layer 1')
        print()
        

        ##################################################################
        # Calculation of Hamiltonian Total Energy (Kinetic + Potential Energy)

        const =  100            # At the legend is added "(x 10^-2)"

        u1 = u1_t[:,::-1,:]     # u1ph
        v1 = v1_t[:,::-1,:]     # u1th

        u2 = u2_t[:,::-1,:]     # u2ph
        v2 = v2_t[:,::-1,:]     # u2th

        h1 = h1_t[:,::-1,:]     # Pseudo-height H
        h2 = h2_t[:,::-1,:]     # Pseudo-height H

        b1 = b1_t[:,::-1,:]     # b1
        b2 = b2_t[:,::-1,:]     # b2    
        
        
        # This was checked with Masoud by the paper 
        #--------------------------------- 
            
        h1_tild                 = 0.5*(h1+H1) + (h2+H2)
        
        energy_layer_1          = (h1+H1)*( 0.5*(u1**2 + v1**2) + h1_tild*(b1+B1) )
        e_1                     = energy_layer_1*const

        energy_layer_1_sum      = np.nansum(energy_layer_1, axis=0)

        #---------------------------------
        
        # h2_tild               = 0.5*(h2+H2)
            
        # energy_layer_2        = (h2+H2)*( 0.5*(u2**2 + v2**2) + h2_tild*(b2+B2) ) 
        # e_2                   = energy_layer_2*const
        
        # energy_layer_2_sum    = np.nansum(energy_layer_2, axis=0)
                
        #--------------------------------- 
        
        # energy_layer_all      = energy_layer_1 + energy_layer_2
        # e_all                 = energy_layer_all*const
        
        # energy_layer_all_timsum = np.nansum(energy_layer_all, axis=0) 
        
        #---------------------------------  
        
        # print('min -->', np.min(e_1), np.min(e_2), np.min(e_all))
        # print('max -->', np.max(e_1), np.max(e_2), np.max(e_all))
        

        ##################################################################

        # Meridional mean 
        
        e_layer_1_mean_glb_inlat      = np.nanmean(np.abs(e_1[:, :, :]), axis=(0,1))
        e_layer_1_mean_glb_inlat_hov  = np.nanmean(np.abs(e_1[:, :, :]), axis=1    )
                
        # Zonal mean
        
        e_layer_1_mean_glb_inlon      = np.nanmean(np.abs(e_1[:, :, :]), axis=(0,2))
        e_layer_1_mean_glb_inlon_hov  = np.nanmean(np.abs(e_1[:, :, :]), axis=2    )
                
        #--------------------------------- 
        
        ##################################################################


        
        print()
        print('# Ploting -->  Meridional Mean of Hamiltonian  (Layer 1)  -->  Hovmöller plot  -->  time x lon - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(28.2, 30, 0.1)

        ax2 = fig.add_subplot(gs[1:, 0])
        
        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'orange',  'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lon[:len(e_layer_1_mean_glb_inlat_hov[0,:])], np.swapaxes(e_layer_1_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lon[:len(e_layer_1_mean_glb_inlat_hov[0,:])], np.swapaxes(e_layer_1_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
        ax2.set_yticklabels([u'', u'120\N{DEGREE SIGN}W', u'60\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'60\N{DEGREE SIGN}E', u'120\N{DEGREE SIGN}E', u''][::-1])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[0], -loc[0]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[1], -loc[1]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)

        plt.title(' Aeolus2.0 - Case {} \n Meridional mean of Hamiltonian in Layer 1'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Hamiltonian_Meridional_Mean_Layer_1.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

        if exists(figname): print('done -->', figname)
                
        

        print()
        print('# Ploting -->  Zonal Mean of Hamiltonian  (Layer 1)  -->  Hovmöller plot  -->  time x lat - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(28.2, 30, 0.1)

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'orange', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lat[:len(e_layer_1_mean_glb_inlon_hov[0,:])], np.swapaxes(e_layer_1_mean_glb_inlon_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lat[:len(e_layer_1_mean_glb_inlon_hov[0,:])], np.swapaxes(e_layer_1_mean_glb_inlon_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax2.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[2], loc[2]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[3], loc[3]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.title(' Aeolus2.0 - Case {} \n Zonal mean of Hamiltonian in Layer 1'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Hamiltonian_Zonal_Mean_Layer_1.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.close()

        if exists(figname): print('done -->', figname)
    
        
    # End of plots_hamiltonian_layer_1
    
    


    plots_hamiltonian_layer_2 = 1


    if plots_hamiltonian_layer_2 == 1:
        
        print('')
        print('> Hamiltonian (Kinetic + Potential Energy) in Layer 2')
        print()
        

        ##################################################################
        # Calculation of Hamiltonian Total Energy (Kinetic + Potential Energy)
        
        const =  100            # At the legend is added "(x 10^-2)"

        u1 = u1_t[:,::-1,:]     # u1ph
        v1 = v1_t[:,::-1,:]     # u1th

        u2 = u2_t[:,::-1,:]     # u2ph
        v2 = v2_t[:,::-1,:]     # u2th

        h1 = h1_t[:,::-1,:]     # Pseudo-height H
        h2 = h2_t[:,::-1,:]     # Pseudo-height H

        b1 = b1_t[:,::-1,:]     # b1
        b2 = b2_t[:,::-1,:]     # b2    
        
        
        # This was checked with Masoud by the paper 
        #--------------------------------- 
            
        # h1_tild               = 0.5*(h1+H1) + (h2+H2)
        
        # energy_layer_1        = (h1+H1)*( 0.5*(u1**2 + v1**2) + h1_tild*(b1+B1) )
        # e_1                   = energy_layer_1*const

        # energy_layer_1_sum    = np.nansum(energy_layer_1, axis=0)

        #---------------------------------
        
        h2_tild                 = 0.5*(h2+H2)
            
        energy_layer_2          = (h2+H2)*( 0.5*(u2**2 + v2**2) + h2_tild*(b2+B2) ) 
        e_2                     = energy_layer_2*const
        
        energy_layer_2_sum      = np.nansum(energy_layer_2, axis=0)
                
        #--------------------------------- 
        
        # energy_layer_all      = energy_layer_1 + energy_layer_2
        # e_all                 = energy_layer_all*const
        
        # energy_layer_all_timsum = np.nansum(energy_layer_all, axis=0) 
        
        #---------------------------------  
        
        # print('min -->', np.min(e_1), np.min(e_2), np.min(e_all))
        # print('max -->', np.max(e_1), np.max(e_2), np.max(e_all))
        

        ##################################################################

        # Meridional mean 
                
        e_layer_2_mean_glb_inlat      = np.nanmean(np.abs(e_2[:, :, :]), axis=(0,1))
        e_layer_2_mean_glb_inlat_hov  = np.nanmean(np.abs(e_2[:, :, :]), axis=1    )
        
        # Zonal mean
                
        e_layer_2_mean_glb_inlon      = np.nanmean(np.abs(e_2[:, :, :]), axis=(0,2))
        e_layer_2_mean_glb_inlon_hov  = np.nanmean(np.abs(e_2[:, :, :]), axis=2    )
        
        #--------------------------------- 
        
        ##################################################################


    
        print()
        print('# Ploting -->  Meridional Mean of Hamiltonian  (Layer 2)  -->  Hovmöller plot  -->  time x lon - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(22.9, 23.8, 0.1)

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'orange',  'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lon[:len(e_layer_2_mean_glb_inlat_hov[0,:])], np.swapaxes(e_layer_2_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lon[:len(e_layer_2_mean_glb_inlat_hov[0,:])], np.swapaxes(e_layer_2_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
        ax2.set_yticklabels([u'', u'120\N{DEGREE SIGN}W', u'60\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'60\N{DEGREE SIGN}E', u'120\N{DEGREE SIGN}E', u''][::-1])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[0], -loc[0]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[1], -loc[1]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)

        plt.title(' Aeolus2.0 - Case {} \n Meridional mean of Hamiltonian in Layer 2'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Hamiltonian_Meridional_Mean_Layer_2.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf()

        if exists(figname): print('done -->', figname)
                
        

        print()
        print('# Ploting -->  Zonal Mean of Hamiltonian  (Layer 2)  -->  Hovmöller plot  -->  time x lat - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs = np.arange(22.9, 23.8, 0.1)

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['white', 'mistyrose', 'salmon', 'orange', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lat[:len(e_layer_2_mean_glb_inlon_hov[0,:])], np.swapaxes(e_layer_2_mean_glb_inlon_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lat[:len(e_layer_2_mean_glb_inlon_hov[0,:])], np.swapaxes(e_layer_2_mean_glb_inlon_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='max')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax2.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[2], loc[2]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[3], loc[3]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.title(' Aeolus2.0 - Case {} \n Zonal mean of Hamiltonian in Layer 2'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_Hamiltonian_Zonal_Mean_Layer_2.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.close()

        if exists(figname): print('done -->', figname)

        
    # End of plots_hamiltonian_layer_2

    
    
        
    plots_wind_layer_1 = 1


    if plots_wind_layer_1 == 1:
        
        print('')
        print('> Wind U and V in Layer 1')
        print()
        

        ##################################################################
        # Wind
        
        const =  100            # At the legend is added "(x 10^-2)"

        u1 = u1_t[:,:,:]        # u1ph
        v1 = v1_t[:,:,:]        # u1th

        #--------------------------------- 
        # Meridional mean
        
        v1_mean_glb_inlat       = np.nanmean(v1[:, :, :], axis=(0,1))*const
        v1_mean_glb_inlat_hov   = np.nanmean(v1[:, :, :], axis=1    )*const

        #--------------------------------- 
        # Zonal mean
        
        u1_mean_glb_inlon       = np.nanmean(u1[:, :, :], axis=(0,2))*const
        u1_mean_glb_inlon_hov   = np.nanmean(u1[:, :, :], axis=2    )*const
        
        ##################################################################


        
        print()
        print('# Ploting -->  Meridional Mean of V Wind (Layer 1)  -->  Hovmöller plot  -->  time x lon - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs1 = np.arange(-0.6, -0.02, 0.04)
        clevs2 = np.arange(0.02, 0.6+0.01, 0.04)
        clevs  = np.concatenate([clevs1, clevs2])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['lime', 'blue', 'mediumblue', 'cornflowerblue', 'lightsteelblue', 'white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lon[:len(v1_mean_glb_inlat_hov[0,:])], np.swapaxes(v1_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lon[:len(v1_mean_glb_inlat_hov[0,:])], np.swapaxes(v1_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='both')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
        ax2.set_yticklabels([u'', u'120\N{DEGREE SIGN}W', u'60\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'60\N{DEGREE SIGN}E', u'120\N{DEGREE SIGN}E', u''][::-1])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[0], -loc[0]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[1], -loc[1]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)

        plt.title(' Aeolus2.0 - Case {} \n Meridional mean of V wind in Layer 1'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_V_Meridional_Mean_Layer_1.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

        if exists(figname): print('done -->', figname)
                
        

        print()
        print('# Ploting -->  Zonal Mean of U Wind (Layer 1)  -->  Hovmöller plot  -->  time x lat - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs1 = np.arange(-0.6, -0.02, 0.04)
        clevs2 = np.arange(0.02, 0.6+0.01, 0.04)
        clevs  = np.concatenate([clevs1, clevs2])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['lime', 'blue', 'mediumblue', 'cornflowerblue', 'lightsteelblue', 'white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lat[:len(u1_mean_glb_inlon_hov[0,:])], np.swapaxes(u1_mean_glb_inlon_hov, 0,1), clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lat[:len(u1_mean_glb_inlon_hov[0,:])], np.swapaxes(u1_mean_glb_inlon_hov, 0,1), clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='both')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax2.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[2], loc[2]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[3], loc[3]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.title(' Aeolus2.0 - Case {} \n Zonal mean of U wind in Layer 1'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_U_Zonal_Mean_Layer_1.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.close()

        if exists(figname): print('done -->', figname)
        
        
    # End of plots_wind_layer_1

    
    
       
    plots_wind_layer_2 = 1


    if plots_wind_layer_2 == 1:
        
        print('')
        print('> Wind U and V in Layer 2')
        print()
        

        ##################################################################
        # Wind

        const =  100            # At the legend is added "(x 10^-2)"
        
        u2 = u2_t[:,:,:]        # u2ph
        v2 = v2_t[:,:,:]        # u2th

        #--------------------------------- 
        # Meridional mean
        
        v2_mean_glb_inlat       = np.nanmean(v2[:, :, :], axis=(0,1))*const
        v2_mean_glb_inlat_hov   = np.nanmean(v2[:, :, :], axis=1    )*const

        #--------------------------------- 
        # Zonal mean
        
        u2_mean_glb_inlon       = np.nanmean(u2[:, :, :], axis=(0,2))*const
        u2_mean_glb_inlon_hov   = np.nanmean(u2[:, :, :], axis=2    )*const
        
        ##################################################################


        
        print()
        print('# Ploting -->  Meridional Mean of V Wind (Layer 2)  -->  Hovmöller plot  -->  time x lon - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs1 = np.arange(-0.6, -0.02, 0.04)
        clevs2 = np.arange(0.02, 0.6+0.01, 0.04)
        clevs  = np.concatenate([clevs1, clevs2])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['lime', 'blue', 'mediumblue', 'cornflowerblue', 'lightsteelblue', 'white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lon[:len(v2_mean_glb_inlat_hov[0,:])], np.swapaxes(v2_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lon[:len(v2_mean_glb_inlat_hov[0,:])], np.swapaxes(v2_mean_glb_inlat_hov, 0,1)[::-1, :], clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='both')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-180, -120, -60, 0, 60, 120, 180])
        ax2.set_yticklabels([u'', u'120\N{DEGREE SIGN}W', u'60\N{DEGREE SIGN}W', u'0\N{DEGREE SIGN}', u'60\N{DEGREE SIGN}E', u'120\N{DEGREE SIGN}E', u''][::-1])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[0], -loc[0]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [-loc[1], -loc[1]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-180, 180], '--', color='black', alpha=0.8, linewidth=2)

        plt.title(' Aeolus2.0 - Case {} \n Meridional mean of V wind in Layer 2'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_V_Meridional_Mean_Layer_2.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.clf() ; plt.close()

        if exists(figname): print('done -->', figname)
                
        

        print()
        print('# Ploting -->  Zonal Mean of U Wind (Layer 2)  -->  Hovmöller plot  -->  time x lat - Horizontal')
        print()
                
        # Create the Hovmöller plot
        projection = ccrs.PlateCarree()

        fig = plt.figure(figsize=(8, 5), constrained_layout=False)

        # define space for the two plots
        gs = fig.add_gridspec(nrows=4, ncols=1, hspace=0)
        
        # create lower Hovmoeller diagram
        clevs1 = np.arange(-0.6, -0.02, 0.04)
        clevs2 = np.arange(0.02, 0.6+0.01, 0.04)
        clevs  = np.concatenate([clevs1, clevs2])

        ax2 = fig.add_subplot(gs[1:, 0])

        cmap_new = clr.LinearSegmentedColormap.from_list('bump_colors', ['lime', 'blue', 'mediumblue', 'cornflowerblue', 'lightsteelblue', 'white', 'mistyrose', 'salmon', 'firebrick', 'red', 'yellow'], N=256)

        cf = ax2.contourf(dtime_hov, lat[:len(u2_mean_glb_inlon_hov[0,:])], np.swapaxes(u2_mean_glb_inlon_hov, 0,1), clevs, cmap=cmap_new, extend='both')

        cs = ax2.contour (dtime_hov, lat[:len(u2_mean_glb_inlon_hov[0,:])], np.swapaxes(u2_mean_glb_inlon_hov, 0,1), clevs, colors='k', linewidths=0.25)

        cbar = plt.colorbar(cf, orientation='vertical', pad=0.01, shrink=0.8, aspect=40, extendrect=False, extend='both')
        cbar.set_label('(x $10^{-2}$)')

        ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax2.set_yticklabels([u'90\N{DEGREE SIGN}S', u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S', u'0\N{DEGREE SIGN}', u'30\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N', u'90\N{DEGREE SIGN}N'])
        ax2.set_xticks( dtime_hov[4::16])
        ax2.set_xticklabels(times[4::16], rotation=90)
        
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[2], loc[2]], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[0  ], dtime_hov[-1 ]], [loc[3], loc[3]], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.plot([dtime_hov[21 ], dtime_hov[21 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[26 ], dtime_hov[26 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[35 ], dtime_hov[35 ]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[111], dtime_hov[111]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        plt.plot([dtime_hov[199], dtime_hov[199]], [-90, 90], '--', color='black', alpha=0.8, linewidth=2)
        
        plt.title(' Aeolus2.0 - Case {} \n Zonal mean of U wind in Layer 2'.format(exp.replace('_',' ')), fontsize=15, pad=15)
        
        # Saving #########################################################

        figname = '{}/Aeolus2.0_Case_{}_Hovmoller_U_Zonal_Mean_Layer_2.png'.format(dir_figs, exp)

        fig.savefig(figname, format='png', dpi=200, bbox_inches='tight') ; plt.close()

        if exists(figname): print('done -->', figname)
        
        
    # End of plots_wind_layer_2
    
    
    
    
    
    
    



print('')
print('')
print('')
print('')
print('')
print('The End')
print('')
print('')


