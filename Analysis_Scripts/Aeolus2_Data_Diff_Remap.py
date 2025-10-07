################################################################################
#
# Author: Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Date:   01.08.2025
# Type:   Python3 + CDO
#
# Description:
# Script to compute difference between Bump and Control cases and regrid to ~1 degree (360x180).
#
# Note: 
# Given the amount of data, the simulations' output were regridded from ~0.5 degrees to ~1 degree using CDO "remapbil,r360x180".
# The data was checked after the interpolation, where, besides small differences, it has no significant deviation from the original.
#
################################################################################

import os
import sys
import argparse
import warnings
import humanize
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('--overwrite', help='just do --overwrite', action='store_true')
args = parser.parse_args()

overwrite = args.overwrite

def   size(f): return humanize.naturalsize(os.path.getsize(f))
def exists(f): return os.path.exists(f)
	
	
print('')
print('Aeolus2.0 Analysis Preparation - Diff and Remap')
print('')


def cases_diff(fin_exp, fin_control, fout, overwrite=0):

    """
    Function to calculate nc difference.
    """

    print()
    print()
    
    for var in ['u1ph', 'u1th', 'u2ph', 'u2th', 'h1', 'h2', 'b1', 'b2', 'CC1', 'w2']:

        fout_v = fout.replace('.nc','_{}.nc'.format(var))

        if exists(fout_v) and overwrite:
            print('exists -->', fout_v, size(fout_v), '--> will be overwritten')
            
        if exists(fout_v) and not overwrite:
            print('exists -->', fout_v, size(fout_v), '--> will not be overwritten')
            continue

        p_0 = 'sellonlatbox,-180,180,-90,90'
        p_1 = '-chname,lon,longitude,lat,latitude'
        p_2 = '-selname,{}'.format(var)
        p_3 = '-seldate,1980-06-01T06:00:00,1980-08-05T18:00:00' 
        p_4 = '-sub'
        
        cmd = 'cdo -O -f nc4 -z zip {0} {1} {2} {3} {4}  {5} {6} {7}'.format(p_0, p_1, p_2, p_3, p_4, fin_exp, fin_control, fout_v)
    
        print(cmd)

        os.system(cmd)
        
        print()
        
        if exists(fout_v):
            print('done -->', fout_v, size(fout_v))
        


def cases_remap(fin, fout, overwrite=0):

    """
    Function to interpolate nc.
    """

    print()
    print()

    for var in ['u1ph', 'u1th', 'u2ph', 'u2th', 'h1', 'h2', 'b1', 'b2', 'CC1', 'w2']:

        fin_v  =  fin.replace('.nc','_{}.nc'.format(var))
        fout_v = fout.replace('_remapbil_1dg.nc','_{}_remapbil_1dg.nc'.format(var))
                
        if exists(fout_v) and overwrite:
            print('exists -->', fout_v, size(fout_v), '--> will be overwritten')
            
        if exists(fout_v) and not overwrite:
            print('exists -->', fout_v, size(fout_v), '--> will not be overwritten')
            continue

        p_0 = 'sellonlatbox,-180,180,-90,90'
        p_1 = '-chname,lon,longitude,lat,latitude'
        p_2 = '-remapbil,r360x180'
        
        cmd = 'cdo -O -f nc4 -z zip {0} {1} {2} {3} {4}'.format(p_0, p_1, p_2, fin_v, fout_v)
        
        print(cmd)

        os.system(cmd)
        
        print()
        
        if exists(fout_v):
            print('done -->', fout_v, size(fout_v))



for exp in [
'../Data/Aeolus2.0_Output_Dry_Barotropic_Weak/Aeolus2.0_Output_Dry_Barotropic_Weak.nc',
'../Data/Aeolus2.0_Output_Dry_Barotropic_Strong/Aeolus2.0_Output_Dry_Barotropic_Strong.nc',
'../Data/Aeolus2.0_Output_Dry_Baroclinic_Weak/Aeolus2.0_Output_Dry_Baroclinic_Weak.nc',
'../Data/Aeolus2.0_Output_Dry_Baroclinic_Strong/Aeolus2.0_Output_Dry_Baroclinic_Strong.nc',
'../Data/Aeolus2.0_Output_MC_Barotropic_Weak/Aeolus2.0_Output_MC_Barotropic_Weak.nc',
'../Data/Aeolus2.0_Output_MC_Barotropic_Strong/Aeolus2.0_Output_MC_Barotropic_Strong.nc',
'../Data/Aeolus2.0_Output_MC_Baroclinic_Weak/Aeolus2.0_Output_MC_Baroclinic_Weak.nc',
'../Data/Aeolus2.0_Output_MC_Baroclinic_Strong/Aeolus2.0_Output_MC_Baroclinic_Strong.nc']:

    if '_Dry_' in exp:
          
        control = '../Data/Aeolus2.0_Output_Dry_Control/Aeolus2.0_Output_Dry_Control.nc'
        diff    = '{}_minus_Dry_Control.nc'.format(exp[:-3])
        
    if '_MC_'  in exp:
          
        control = '../Data/Aeolus2.0_Output_MC_Control/Aeolus2.0_Output_MC_Control.nc'
        diff    = '{}_minus_MC_Control.nc'.format(exp[:-3])
        
    # Difference
    
    cases_diff(exp, control, diff, overwrite=overwrite)
    
    
    # Remap
    
    diff_remap = '{}_remapbil_1dg.nc'.format(diff[:-3])
    
    cases_remap(diff, diff_remap, overwrite=overwrite)
    












