################################################################################
#
# Author:    Sullyandro Guimaraes (sullyandro@pik-potsdam.de)
# Date:      13.12.2025
# Type:      Python3 
#
# Description:
# Script to download ERA5 data.
#
# This is a adaptation of the API request available on the page:
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels
#
################################################################################

import os
import cdsapi  	
import humanize

def   size(f): return humanize.naturalsize(os.path.getsize(f))
def exists(f): return os.path.exists(f)

if not exists('../Data/ERA5/'): os.makedirs('../Data/ERA5/')

target = '../Data/ERA5/ERA5_6hours_1980-06-01-00_1980-06-30-18.nc'


if exists(target): 
	print('exists -->', target, size(target)) ; exit()
			
dataset = 'reanalysis-era5-pressure-levels'
request = {
    'product_type':   ['reanalysis'],
    'variable':       ['temperature', 'u_component_of_wind', 'v_component_of_wind'],
    'year':           ['1980'],
    'month':          ['06'  ],
    'day':            ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30'],
    'time':           ['00:00', '06:00', '12:00', '18:00'],
    'pressure_level': ['250', '450', '650', '800', '900', '1000'],
    'data_format':    'netcdf',
    'download_format':'unarchived'
}

client = cdsapi.Client()
client.retrieve(dataset, request, target)

print()
if exists(target):
    print('done -->', target, size(target))
    
