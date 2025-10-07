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


print('')
print('Aeolus2.0 Analysis Preparation - Diff and Remap')
print('')


def cases_diff(fin_exp, fin_control, fout, overwrite=0):

    """
    Function to calculate nc difference.
    """

    if os.path.exists(fout) and overwrite:
        print('exists -->', fout, humanize.naturalsize(os.path.getsize(fout)), '--> will be overwritten')
        
    if os.path.exists(fout) and not overwrite:
        print('exists -->', fout, humanize.naturalsize(os.path.getsize(fout)), '--> will not be overwritten')
        return

    print()
    print()
    
    p_0 = 'sellonlatbox,-180,180,-90,90'
    p_1 = '-selname,u1ph,u1th,u2ph,u2th,h1,h2,b1,b2,CC1,w2'
    p_2 = '-seldate,1980-06-01T06:00:00,1980-08-05T18:00:00' 
    p_3 = '-sub'
    
    cmd = 'cdo -O -f nc4 -z zip chname,lon,longitude,lat,latitude {0} {1} {2} {3}  {4} {5} {6}'.format(p_0, p_1, p_2, p_3, fin_exp, fin_control, fout)
    
    print(cmd)

    os.system(cmd)
    
    print()
    
    if os.path.exists(fout):
        print('done -->', fout, humanize.naturalsize(os.path.getsize(fout)))


def cases_remap(fin, fout, overwrite=0):

    """
    Function to interpolate nc.
    """

    if os.path.exists(fout) and overwrite:
        print('exists -->', fout, humanize.naturalsize(os.path.getsize(fout)), '--> will be overwritten')
        
    if os.path.exists(fout) and not overwrite:
        print('exists -->', fout, humanize.naturalsize(os.path.getsize(fout)), '--> will not be overwritten')
        return

    print()
    print()

    p_0 = 'sellonlatbox,-180,180,-90,90'
    p_1 = '-chname,lon,longitude,lat,latitude'
    p_2 = '-remapbil,r360x180'
    
    cmd = 'cdo -O -f nc4 -z zip {0} {1} {2} {3} {4}'.format(p_0, p_1, p_2, fin, fout)
    
    print(cmd)

    os.system(cmd)
    
    print()
    
    if os.path.exists(fout):
        print('done -->', fout, humanize.naturalsize(os.path.getsize(fout)))



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
	












