################################################################################
#
# Author:       Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Colaborators: Masoud Rostami
# Date:         01.08.2025
# Type:         Python3
#
# Description:
# Script to plot maps for Aeolus2 analysis of Bump evolution.
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
print('Aeolus2.0 Analysis Preparation - Plot Maps')


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
    print('# Maps production')
    print()
    



    if remap == 0: dir_figs = '../Figures/Plot_Maps/Aeolus2.0_{}'.format(exp)
    if remap == 1: dir_figs = '../Figures_remapbil_1dg/Plot_Maps/Aeolus2.0_{}'.format(exp)
    
    if not exists(dir_figs): os.makedirs(dir_figs)




    # Defining plot coordenates
    lats_plot = lat[:]
    lons_plot = np.append(lon, [180]) # add cyclic points for plot continuity
    
    # make 2d grid of lons x lats
    lons, lats = np.meshgrid(lons_plot, lats_plot)
    
    # define parallels and meridians to draw
    parallels = np.arange(-80,  90, 20)
    meridians = np.arange(  0, 360, 20)
    
    
    
        
    maps_b_v = 1


    if maps_b_v == 1:
        
        print('')
        print('> Buoyancy and Wind Vectors')
        print()
        
        for i in [21, 26, 35, 111, 199]:

            print('# Ploting -->', i, dtime[i])
            
            ##################################################################
            
            const =  100            # to improve visibility of wind arrows --> At the legend is added "(x 10^-2)"

            u1 = u1_t[i,:,:]*const  # u1ph
            v1 = v1_t[i,:,:]*const  # u1th

            u2 = u2_t[i,:,:]*const  # u2ph
            v2 = v2_t[i,:,:]*const  # u2th
            
            b1 = b1_t[i,:,:]*const  # b1
            b2 = b2_t[i,:,:]*const  # b2


            # applying filter to make the plot soft
            
            u1f = uniform_filter(u1[:,:], size=5, mode='constant')
            v1f = uniform_filter(v1[:,:], size=5, mode='constant')
            
            u2f = uniform_filter(u2[:,:], size=5, mode='constant')
            v2f = uniform_filter(v2[:,:], size=5, mode='constant')


            # add cyclic points manually

            b1n = np.append( b1[:,:],  b1[:,0:1], axis=1)
            u1n = np.append(u1f[:,:], u1f[:,0:1], axis=1)
            v1n = np.append(v1f[:,:], v1f[:,0:1], axis=1)
            
            b2n = np.append( b2[:,:],  b2[:,0:1], axis=1)
            u2n = np.append(u2f[:,:], u2f[:,0:1], axis=1)
            v2n = np.append(v2f[:,:], v2f[:,0:1], axis=1)

            ##################################################################

                                
            fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))


            # Plotting Layer 1 ###############################################
            
            plt.subplot(1, 2, 1)
                        
            # make orthographic basemap
            m = Basemap(resolution='c', projection='ortho', lat_0=lat_0, lon_0=lon_0)

            # set desired contour levels
            clevs = np.linspace(-7.1, 7.1, 10)

            # compute native x,y coordinates of grid
            x, y = m(lons, lats)

            # plot contours
            CS1 = m.contour (x, y, b1n, clevs, linewidths=0.5, colors='k')
            CS2 = m.contourf(x, y, b1n, clevs, cmap=plt.cm.RdBu_r, extend='both')

            # transform vectors to projection grid
            uproj, vproj, xx, yy = m.transform_vector(u1n, v1n, lons_plot, lats_plot, 31, 31, returnxy=True, masked=True)
    
            # now plot
            Q = m.quiver(xx[:], yy[:], uproj[:], vproj[:], scale=const)

            # draw parallels and meridians
            m.drawcoastlines(linewidth=1.5) 
            m.drawparallels(parallels)      
            m.drawmeridians(meridians)      

            # add colorbar
            kwargs = {'format': '%.1f'}
            cb = m.colorbar(CS2, 'bottom', size='5%', pad='2%', **kwargs)
            cb.set_label('(x $10^{-2}$)')

            # Add circle
            xc1, yc1 = m(lon_0, lat_0)
            xc2, yc2 = m(lon_0, lat_0+12) 
            circle1  = plt.Circle((xc1, yc1), yc2-yc1, ls='--', lw=5, color='green', fill=False, alpha=0.5)
            ax1.add_patch(circle1)

            plt.title(' Aeolus2.0 - Case {} \n Buoyancy and Wind Vectors in Layer 1 \n {} '.format(exp.replace('_',' '), str(dtime[i])), fontsize=15)

            
            # Plotting Layer 2 ###############################################
            
            plt.subplot(1, 2, 2)
            
            # make orthographic basemap
            m = Basemap(resolution='c', projection='ortho', lat_0=lat_0, lon_0=lon_0)

            # set desired contour levels
            clevs = np.linspace(-7.1, 7.1, 10)

            # compute native x,y coordinates of grid
            x, y = m(lons, lats)

            # plot contours
            CS1 = m.contour (x, y, b2n, clevs, linewidths=0.5, colors='k')
            CS2 = m.contourf(x, y, b2n, clevs, cmap=plt.cm.RdBu_r, extend='both')
            
            # transform vectors to projection grid
            uproj, vproj, xx, yy = m.transform_vector(u2n, v2n, lons_plot, lats_plot, 31, 31, returnxy=True, masked=True)
    
            # now plot
            Q = m.quiver(xx[:], yy[:], uproj[:], vproj[:], scale=const)

            # draw parallels and meridians
            m.drawcoastlines(linewidth=1.5)
            m.drawparallels(parallels)      
            m.drawmeridians(meridians)      

            # add colorbar
            kwargs = {'format': '%.1f'}
            cb = m.colorbar(CS2, 'bottom', size='5%', pad='2%', **kwargs)
            cb.set_label('(x $10^{-2}$)')

            # Add circle
            xc1, yc1 = m(lon_0, lat_0)
            xc2, yc2 = m(lon_0, lat_0+12) 
            circle2  = plt.Circle((xc1, yc1), yc2-yc1, ls='--', lw=5, color='green', fill=False, alpha=0.5)
            ax2.add_patch(circle2)

            plt.title(' Aeolus2.0 - Case {} \n Buoyancy and Wind Vectors in Layer 2 \n {} '.format(exp.replace('_',' '), str(dtime[i])), fontsize=15)

            
            # Saving #########################################################

            figname = '{}/Aeolus2.0_Case_{}_Buoyancy_and_Wind_Vectors_Step_{:03d}.png'.format(dir_figs, exp, i)

            fig1.savefig(figname, format='png', dpi=80, bbox_inches='tight') ; plt.clf()

            if exists(figname): print('done -->', figname)
            
                        
    # End of maps_b_v
    
    

        
    maps_clwc_w2_v = 1


    if maps_clwc_w2_v == 1 and 'Dry' not in exp:
        
        print('')
        print('> CLWC, W2 and Wind Vectors')
        print()
                
        for i in [21, 26, 35, 111, 199]:

            print('# Ploting -->', i, dtime[i])
            
            ##################################################################
            
            const =  100               # to improve visibility of wind arrows --> At the legend is added "(x 10^-2)"

            u1   = u1_t[i,:,:]*const   # u1ph
            v1   = v1_t[i,:,:]*const   # u1th

            u2   = u2_t[i,:,:]*const   # u2ph
            v2   = v2_t[i,:,:]*const   # u2th
            
            clwc = cc1_t[i,:,:]*const  # Condensed liquid water content (CLWC)
            w2   =  w2_t[i,:,:]*const  # Bulk of precipitable water


            # applying filter to make the plot soft
            
            u1f = uniform_filter(u1[:,:], size=5, mode='constant')
            v1f = uniform_filter(v1[:,:], size=5, mode='constant')
            
            u2f = uniform_filter(u2[:,:], size=5, mode='constant')
            v2f = uniform_filter(v2[:,:], size=5, mode='constant')


            # add cyclic points manually

            clwcn = np.append(clwc[:,:], clwc[:,0:1], axis=1)
            u1n   = np.append(u1f[:,:],   u1f[:,0:1], axis=1)
            v1n   = np.append(v1f[:,:],   v1f[:,0:1], axis=1)
            
            w2n   = np.append( w2[:,:],    w2[:,0:1], axis=1)
            u2n   = np.append(u2f[:,:],   u2f[:,0:1], axis=1)
            v2n   = np.append(v2f[:,:],   v2f[:,0:1], axis=1)

            ##################################################################

                                
            fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))


            # Plotting Layer 1 ###############################################
            
            plt.subplot(1, 2, 1)
                        
            # make orthographic basemap
            m = Basemap(resolution='c', projection='ortho', lat_0=lat_0, lon_0=lon_0)

            # set desired contour levels
            clevs = np.linspace(-0.17, 0.17, 20)

            # compute native x,y coordinates of grid
            x, y = m(lons, lats)

            # plot contours
            CS1 = m.contour (x, y, clwcn, clevs, linewidths=0.5, colors='k')
            CS2 = m.contourf(x, y, clwcn, clevs, cmap=plt.cm.RdBu, extend='both')

            # transform vectors to projection grid
            uproj, vproj, xx, yy = m.transform_vector(u1n, v1n, lons_plot, lats_plot, 31, 31, returnxy=True, masked=True)
    
            # now plot
            Q = m.quiver(xx[:], yy[:], uproj[:], vproj[:], scale=const)

            # draw parallels and meridians
            m.drawcoastlines(linewidth=1.5)
            m.drawparallels(parallels)      
            m.drawmeridians(meridians)     

            # add colorbar
            kwargs = {'format': '%.2f'}
            cb = m.colorbar(CS2, 'bottom', size='5%', pad='2%', **kwargs)
            cb.set_label('(x $10^{-2}$)')

            # Add circle
            xc1, yc1 = m(lon_0, lat_0)
            xc2, yc2 = m(lon_0, lat_0+12)
            circle1  = plt.Circle((xc1, yc1), yc2-yc1, ls='--', lw=5, color='green', fill=False, alpha=0.5)
            ax1.add_patch(circle1)

            plt.title(' Aeolus2.0 - Case {} \n CLWC and Wind Vectors in Layer 1 \n {} '.format(exp.replace('_',' '), str(dtime[i])), fontsize=15)

            
            # Plotting Layer 2 ###############################################
            
            plt.subplot(1, 2, 2)
            
            # make orthographic basemap
            m = Basemap(resolution='c', projection='ortho', lat_0=lat_0, lon_0=lon_0)

            # set desired contour levels
            clevs = np.linspace(-0.17, 0.17, 20)

            # compute native x,y coordinates of grid
            x, y = m(lons, lats)

            # plot contours
            CS1 = m.contour (x, y, w2n, clevs, linewidths=0.5, colors='k')
            CS2 = m.contourf(x, y, w2n, clevs, cmap=plt.cm.RdBu, extend='both')
            
            # transform vectors to projection grid
            uproj, vproj, xx, yy = m.transform_vector(u2n, v2n, lons_plot, lats_plot, 31, 31, returnxy=True, masked=True)
    
            # now plot
            Q = m.quiver(xx[:], yy[:], uproj[:], vproj[:], scale=const)

            # draw parallels and meridians
            m.drawcoastlines(linewidth=1.5) #, color='gray')
            m.drawparallels(parallels)      #, color='gray')
            m.drawmeridians(meridians)      #, color='gray')

            # add colorbar
            kwargs = {'format': '%.2f'}
            cb = m.colorbar(CS2, 'bottom', size='5%', pad='2%', **kwargs)
            cb.set_label('(x $10^{-2}$)')

            # Add circle
            xc1, yc1 = m(lon_0, lat_0)
            xc2, yc2 = m(lon_0, lat_0+12) 
            circle2  = plt.Circle((xc1, yc1), yc2-yc1, ls='--', lw=5, color='green', fill=False, alpha=0.5)
            ax2.add_patch(circle2)

            plt.title(' Aeolus2.0 - Case {} \n W2 and Wind Vectors in Layer 2 \n {} '.format(exp.replace('_',' '), str(dtime[i])), fontsize=15)

            
            # Saving #########################################################

            figname = '{}/Aeolus2.0_Case_{}_CLWC_W2_and_Wind_Vectors_Step_{:03d}.png'.format(dir_figs, exp, i)

            fig1.savefig(figname, format='png', dpi=80, bbox_inches='tight') ; plt.clf()

            if exists(figname): print('done -->', figname)
            
                    
    # End of maps_clwc_w2_v
    
    
    
        
    maps_hamiltonian = 1


    if maps_hamiltonian == 1:
        
        print('')
        print('> Hamiltonian (Kinetic + Potential Energy)')
        print()
        
        for i in [21, 26, 35, 111, 199]:

            print('# Ploting -->', i, dtime[i])


            ##################################################################
            # Calculation of Hamiltonian Total Energy (Kinetic + Potential Energy)

            const =  100        # to improve visibility of wind arrows --> At the legend is added "(x 10^-2)"

            u1 = u1_t[i,:,:]    # u1ph
            v1 = v1_t[i,:,:]    # u1th

            u2 = u2_t[i,:,:]    # u2ph
            v2 = v2_t[i,:,:]    # u2th

            h1 = h1_t[i,:,:]    # Pseudo-height H
            h2 = h2_t[i,:,:]    # Pseudo-height H

            b1 = b1_t[i,:,:]    # b1
            b2 = b2_t[i,:,:]    # b2    
            
            
            # This was checked with Masoud by the paper 
            #--------------------------------- 
                
            h1_tild                 = 0.5*(h1+H1) + (h2+H2)
            
            energy_layer_1          = (h1+H1)*( 0.5*(u1**2 + v1**2) + h1_tild*(b1+B1) )
            e_1                     = energy_layer_1*const

            energy_layer_1_sum      = np.nansum(energy_layer_1, axis=0)

            #---------------------------------
            
            h2_tild                 = 0.5*(h2+H2)
                
            energy_layer_2          = (h2+H2)*( 0.5*(u2**2 + v2**2) + h2_tild*(b2+B2) ) 
            e_2                     = energy_layer_2*const
            
            energy_layer_2_sum      = np.nansum(energy_layer_2, axis=0)
                    
            #--------------------------------- 
            
            energy_layer_all        = energy_layer_1 + energy_layer_2
            e_all                   = energy_layer_all*const
            
            energy_layer_all_timsum = np.nansum(energy_layer_all, axis=0) 
            
            #---------------------------------  
            
            # print('min -->', np.min(e_1), np.min(e_2), np.min(e_all))
            # print('max -->', np.max(e_1), np.max(e_2), np.max(e_all))
            

            # add cyclic points manually

            e_1n = np.append( e_1[:,:],  e_1[:,0:1], axis=1)
            e_2n = np.append( e_2[:,:],  e_2[:,0:1], axis=1)

            ##################################################################

                                
            fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))


            # Plotting Layer 1 ###############################################
            
            plt.subplot(1, 2, 1)
                        
            # make orthographic basemap
            m = Basemap(resolution='c', projection='ortho', lat_0=lat_0, lon_0=lon_0)

            # set desired contour levels
            clevs = np.linspace(-50, 50, 20)

            # compute native x,y coordinates of grid
            x, y = m(lons, lats)

            # plot contours
            CS1 = m.contour (x, y, e_1n, clevs, linewidths=0.5, colors='k')
            CS2 = m.contourf(x, y, e_1n, clevs, cmap=plt.cm.RdBu_r, extend='both')

            # draw parallels and meridians
            m.drawcoastlines(linewidth=1.5)
            m.drawparallels(parallels)      
            m.drawmeridians(meridians)      

            # add colorbar
            kwargs = {'format': '%.1f'}
            cb = m.colorbar(CS2, 'bottom', size='5%', pad='2%', **kwargs)
            cb.set_label('(x $10^{-2}$)')

            # Add circle
            xc1, yc1 = m(lon_0, lat_0)
            xc2, yc2 = m(lon_0, lat_0+12) 
            circle1  = plt.Circle((xc1, yc1), yc2-yc1, ls='--', lw=5, color='green', fill=False, alpha=0.5)
            ax1.add_patch(circle1)

            plt.title(' Aeolus2.0 - Case {} \n Hamiltonian (Kinetic + Potential Energy) in Layer 1 \n {} '.format(exp.replace('_',' '), str(dtime[i])), fontsize=15)

            
            # Plotting Layer 2 ###############################################
            
            plt.subplot(1, 2, 2)
            
            # make orthographic basemap
            m = Basemap(resolution='c', projection='ortho', lat_0=lat_0, lon_0=lon_0)

            # set desired contour levels
            clevs = np.linspace(-50, 50, 20)

            # compute native x,y coordinates of grid
            x, y = m(lons, lats)

            # plot contours
            CS1 = m.contour (x, y, e_2n, clevs, linewidths=0.5, colors='k')
            CS2 = m.contourf(x, y, e_2n, clevs, cmap=plt.cm.RdBu_r, extend='both')
            
            # draw parallels and meridians
            m.drawcoastlines(linewidth=1.5)
            m.drawparallels(parallels)      
            m.drawmeridians(meridians)     

            # add colorbar
            kwargs = {'format': '%.1f'}
            cb = m.colorbar(CS2, 'bottom', size='5%', pad='2%', **kwargs)
            cb.set_label('(x $10^{-2}$)')

            # Add circle
            xc1, yc1 = m(lon_0, lat_0)
            xc2, yc2 = m(lon_0, lat_0+12) 
            circle2  = plt.Circle((xc1, yc1), yc2-yc1, ls='--', lw=5, color='green', fill=False, alpha=0.5)
            ax2.add_patch(circle2)

            plt.title(' Aeolus2.0 - Case {} \n Hamiltonian (Kinetic + Potential Energy) in Layer 2 \n {} '.format(exp.replace('_',' '), str(dtime[i])), fontsize=15)

            
            # Saving #########################################################

            figname = '{}/Aeolus2.0_Case_{}_Hamiltonian_Energy_Step_{:03d}.png'.format(dir_figs, exp, i)

            fig1.savefig(figname, format='png', dpi=80, bbox_inches='tight') ; plt.clf()

            if exists(figname): print('done -->', figname)
            
        
    # End of maps_hamiltonian
    
    
    
    

    




print('')
print('')
print('')
print('')
print('')
print('The End')
print('')
print('')


