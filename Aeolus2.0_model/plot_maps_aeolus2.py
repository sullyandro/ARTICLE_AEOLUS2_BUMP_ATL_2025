################################################################################
# Author: Sullyandro Guimaraes (sullyandro@gmail.com)
# Date:   25.01.2024
# Type:   Python3
#
# Description:
# Script to plot maps.
#
################################################################################

"""
netcdf teste_npytonc_output {
dimensions:
	longitude = 768 ;
	latitude = 384 ;
	time = UNLIMITED ; // (123 currently)
variables:
	float longitude(longitude) ;
		longitude:units = "degrees_east" ;
		longitude:axis = "X" ;
		longitude:cartesian_axis = "X" ;
	float latitude(latitude) ;
		latitude:units = "degrees_north" ;
		latitude:axis = "Y" ;
		latitude:cartesian_axis = "Y" ;
	double time(time) ;
		time:units = "days since 0001-01-01 00:00" ;
		time:axis = "T" ;
		time:cartesian_axis = "T" ;
	float h1(time, latitude, longitude) ;
	float u1ph(time, latitude, longitude) ;
	float u1th(time, latitude, longitude) ;
	float h2(time, latitude, longitude) ;
	float RT1(time, latitude, longitude) ;
	float RT2(time, latitude, longitude) ;
	float u2ph(time, latitude, longitude) ;
	float u2th(time, latitude, longitude) ;
	float q1(time, latitude, longitude) ;
	float b1(time, latitude, longitude) ;
	float b2(time, latitude, longitude) ;
	float q2(time, latitude, longitude) ;
	float w1(time, latitude, longitude) ;
	float w2(time, latitude, longitude) ;
	float menthalpy(time, latitude, longitude) ;
	float CC1(time, latitude, longitude) ;
	float Prec1(time, latitude, longitude) ;
	float DD1(time, latitude, longitude) ;
	float DDr(time, latitude, longitude) ;
	float Ev(time, latitude, longitude) ;


yes? us teste_npytonc_output.nc
yes? sh d
     currently SET data sets:
    1> ./teste_npytonc_output.nc  (default)
 name     title                             I         J         K         L         M         N
 H1                                        1:768     1:384     ...       1:123     ...       ...
 U1PH                                      1:768     1:384     ...       1:123     ...       ...
 U1TH                                      1:768     1:384     ...       1:123     ...       ...
 H2                                        1:768     1:384     ...       1:123     ...       ...
 RT1                                       1:768     1:384     ...       1:123     ...       ...
 RT2                                       1:768     1:384     ...       1:123     ...       ...
 U2PH                                      1:768     1:384     ...       1:123     ...       ...
 U2TH                                      1:768     1:384     ...       1:123     ...       ...
 Q1                                        1:768     1:384     ...       1:123     ...       ...
 B1                                        1:768     1:384     ...       1:123     ...       ...
 B2                                        1:768     1:384     ...       1:123     ...       ...
 Q2                                        1:768     1:384     ...       1:123     ...       ...
 W1                                        1:768     1:384     ...       1:123     ...       ...
 W2                                        1:768     1:384     ...       1:123     ...       ...
 MENTHALPY
                                           1:768     1:384     ...       1:123     ...       ...
 CC1                                       1:768     1:384     ...       1:123     ...       ...
 PREC1                                     1:768     1:384     ...       1:123     ...       ...
 DD1                                       1:768     1:384     ...       1:123     ...       ...
 DDR                                       1:768     1:384     ...       1:123     ...       ...
 EV                                        1:768     1:384     ...       1:123     ...       ...

  on grid GRS1 with -1.E+34 for missing data
  X=179.8E(-180.2):179.8E(179.8)  Y=89.9S:89.9N 
             
time range: 01-JAN 00:00 to 05-FEB 05:52

ncout.AddVariable('h10',        'Effective topography height',                 '[gamma_topo*g_real*H0_dim]',               miss_val=None, axes=axes_yx)
ncout.AddVariable('ocean',      'Land-sea mask',                               '1',                                        miss_val=None, axes=axes_yx)
ncout.AddVariable('u1th',       'Meridional velocity layer 1',                 '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
ncout.AddVariable('u1ph',       'Zonal (azimuthal) velocity layer 1',          '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
ncout.AddVariable('h1',         'Pseudo-height layer 1',                       'H',                                        miss_val=None, axes=axes_tyx)
ncout.AddVariable('q1',         'Bulk of Specific humidity at layer 1',        '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
ncout.AddVariable('q2',         'Bulk of Specific humidity at layer 2',        '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
ncout.AddVariable('w1',         'Bulk of Precipitable Water at layer 1',       '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
ncout.AddVariable('w2',         'Bulk of Precipitable Water at layer 2',       '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
ncout.AddVariable('Prec1',      'Bulk of Precipitaion at layer 1',             '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
ncout.AddVariable('b1',         'Buoyancy layer 1',                            '[g*theta/theta_s]',                        miss_val=None, axes=axes_tyx)
ncout.AddVariable('b2',         'Buoyancy layer 2',                            '[g*theta/theta_s]',                        miss_val=None, axes=axes_tyx)
ncout.AddVariable('CC1',        'CLWC (Cloud Liquid Water Content) at layer 1','[Q/T]',                                    miss_val=None, axes=axes_tyx)
ncout.AddVariable('DD1',        'Downdraft to layer 1 (balanced)',             '[Q/T]',                                    miss_val=None, axes=axes_tyx)
ncout.AddVariable('DDr',        'Downdraft to layer 1 (unbalanced)',           '[Q/T]',                                    miss_val=None, axes=axes_tyx)
ncout.AddVariable('Ev',         'Sea surface evaoporation',                    '[Q/T]',                                    miss_val=None, axes=axes_tyx)
ncout.AddVariable('u2th',       'Meridional velocity layer 2',                 '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
ncout.AddVariable('u2ph',       'Zonal (azimuthal) velocity layer 2',          '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
ncout.AddVariable('h2',         'Pseudo-height layer 2',                       'H',                                        miss_val=None, axes=axes_tyx)
ncout.AddVariable('RT1',        'Radiative transfer flux, lower layer',        '[HB]',                                     miss_val=None, axes=axes_tyx)
ncout.AddVariable('RT2',        'Radiative transfer flux, upper layer',        '[HB]',                                     miss_val=None, axes=axes_tyx)
ncout.AddVariable('menthalpy',  'menthalpy = h1+H1-ep1*q1-ep1*Q01',            'unknown units',                            miss_val=None, axes=axes_tyx)
ncout.AddVariable('insolation', 'daily mean insolation',                       'W/m2?',                                    miss_val=None, axes=axes_ty) 


                print('h1,h2                Layers Pseudo-height')
                print('q1                   Bulk of Specific humidity')
                print('m_enthalpy1          Enthalpy')
                print('b1,b2                Buoyancy')
                print('u1,u2                Meridional velocity')
                print('w1                   Bulk of Precipitable Water')
                print('E1                   Sea surface evaoporation')
                print('C1                   CLWC (Cloud Liquid Water Content)')
                print('D1                   Downdraft to layer 1 (balanced)')
                print('nd                   Non dimensional')

"""

import os
import glob
import warnings
import numpy as np
import matplotlib as mpl #; mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from   netCDF4 import Dataset
from   PIL import Image, ImageDraw, ImageFont
import moviepy.editor as mp
plt.style.use('classic')
warnings.filterwarnings("ignore")
# np.set_printoptions(precision=2, suppress=1, linewidth=500)

var_dict = {
'h10'		:{'longname':'Effective topography height',                 'unit':'gamma_topo*g_real*H0_dim',  				},             
'ocean'		:{'longname':'Land-sea mask',                               'unit':'1',                            				},
'u1th'		:{'longname':'Meridional velocity layer 1',                 'unit':'[(g*H)^0.5], [beta*L_d^2] at the equator',	},
'u1ph'		:{'longname':'Zonal (azimuthal) velocity layer 1',          'unit':'[(g*H)^0.5], [beta*L_d^2] at the equator', 	},
'h1'		:{'longname':'Pseudo-height layer 1',                       'unit':'H',                                        	},
'q1'		:{'longname':'Bulk of Specific humidity at layer 1',        'unit':'(L_v.g/(C_p.theta_s))Kg/Kg',             	},
'q2'		:{'longname':'Bulk of Specific humidity at layer 2',        'unit':'(L_v.g/(C_p.theta_s))Kg/Kg',             	},
'w1'		:{'longname':'Bulk of Precipitable Water at layer 1',       'unit':'(L_v.g/(C_p.theta_s))Kg/Kg',             	},
'w2'		:{'longname':'Bulk of Precipitable Water at layer 2',       'unit':'(L_v.g/(C_p.theta_s))Kg/Kg',             	},
'Prec1'		:{'longname':'Bulk of Precipitaion at layer 1',             'unit':'(L_v.g/(C_p.theta_s))Kg/Kg',             	},
'b1'		:{'longname':'Buoyancy layer 1',                            'unit':'g*theta/theta_s',                        	},
'b2'		:{'longname':'Buoyancy layer 2',                            'unit':'g*theta/theta_s',                        	},
'CC1'		:{'longname':'CLWC (Cloud Liquid Water Content) at layer 1','unit':'Q/T',                                    	},
'DD1'		:{'longname':'Downdraft to layer 1 (balanced)',             'unit':'Q/T',                                    	},
'DDr'		:{'longname':'Downdraft to layer 1 (unbalanced)',           'unit':'Q/T',                                    	},
'Ev'		:{'longname':'Sea surface evaoporation',                    'unit':'Q/T',                                    	},
'u2th'		:{'longname':'Meridional velocity layer 2',                 'unit':'[(g*H)^0.5], [beta*L_d^2] at the equator', 	},
'u2ph'		:{'longname':'Zonal (azimuthal) velocity layer 2',          'unit':'[(g*H)^0.5], [beta*L_d^2] at the equator', 	},
'h2'		:{'longname':'Pseudo-height layer 2',                       'unit':'H',                                        	},
'RT1'		:{'longname':'Radiative transfer flux, lower layer',        'unit':'HB',                                     	},
'RT2'		:{'longname':'Radiative transfer flux, upper layer',        'unit':'HB',                                     	},
'menthalpy'	:{'longname':'menthalpy = h1+H1-ep1*q1-ep1*Q01',            'unit':'unknown units',                            	},
'insolation':{'longname':'daily mean insolation',                       'unit':'W/m2?',                                    	} }

print()
print('# Aeolus2 variables:')
print()
print('{:15} {:50} {:50}'.format('Variable', 'Longname', 'Unit'))
print('------------------------------------------------------------------------------------------------------------')
for varname in var_dict:
	print('{:15} {:50} {:50}'.format(varname, var_dict[varname]['longname'], var_dict[varname]['unit']))
print('------------------------------------------------------------------------------------------------------------')
print()

sufix = 'case_june1980_mc0'

print()
print('# Reading NetCDF4 file')

fin = 'output/npytonc_full_output.nc'

dat = Dataset(fin, mode='r')
lon = dat.variables['longitude'][:]
lat = dat.variables['latitude' ][:]
tim = dat.variables['time' ][:]

freq = 4

variables = list(dat.variables.keys())
variables.remove('longitude')
variables.remove('latitude')
variables.remove('time')

print()
print(variables)

print()
print('lats:', lat.shape, '[', lat[0], lat[1], '...', lat[-2], lat[-1], ']')
print()
print('lons:', lon.shape, '[', lon[0], lon[1], '...', lon[-2], lon[-1], ']')


def find_best_levels(sample, number_of_levels=10, percentil_min=1, percentil_max=99):
	"""
	Function to find the best map levels for a variable (sample).
	The variable can be of any shape.
	The idea is to have levels from 5% to 95% every ~10%, for example. 
	"""		
	
	import statsmodels.api as sm

	ecdf = sm.distributions.ECDF(np.reshape(sample, -1))
	
	x = np.linspace(np.min(sample), np.max(sample), 1000)
	y = ecdf(x)
	
	percents = np.linspace(percentil_min, percentil_max, number_of_levels)
		
	index = []
	for i in percents:
		index.append((np.abs(y - i/100)).argmin())
	
	levels = np.unique(x[index])
	
	if len(levels) <= int(number_of_levels/3):
		new_levels = []
		for i in range(len(levels)-1):
			new_levels.append(levels[i])
			new_levels.append(levels[i]+(levels[i+1]-levels[i])*(1/3))
			new_levels.append(levels[i]+(levels[i+1]-levels[i])*(2/3))
		levels = new_levels
		
	if len(levels) <= int(number_of_levels/2):
		new_levels = []
		for i in range(len(levels)-1):
			new_levels.append(levels[i])
			new_levels.append((levels[i]+levels[i+1])/2)
		levels = new_levels
	
	return levels
		
seq_1 = [  'b1',   'q1',    'h1',  
		   'b2',   'q2',    'h2' ]
seq_2 = ['u1ph', 'u1th', 'Prec1', 
		 'u2ph', 'u2th',   'CC1' ]
seq = seq_1 + seq_2

seq = ['Prec1', 'CC1']


print()
print()
print('# Variables loop')
print()
for varname in seq:
	
	longname = var_dict[varname]['longname']
	unit     = var_dict[varname]['unit']
	
	var = dat.variables[varname][:]
	var_timmean     = np.mean(var, axis=0)
	var_timmean_min = np.min(var_timmean)
	var_timmean_max = np.max(var_timmean)
	
	var_timmean_min = np.percentile(var,  3)
	var_timmean_max = np.percentile(var, 97)
	
	limits = max(abs(var_timmean_min), abs(var_timmean_max))
		
	levels = np.linspace(-limits, limits, 10)
	
	print()
	print('Variable -->', varname, var.shape, var_timmean_min, var_timmean_max)
	print()
	
	# levels = find_best_levels(var, number_of_levels=16)
	
	if var_timmean_min == var_timmean_max and varname not in ['Prec1', 'CC1']:
		print()
		print('*** Variable with problem, check it.')
		print()
		continue
	elif var_timmean_min == var_timmean_max and varname in ['Prec1', 'CC1']:
		var_timmean_min = -1
		var_timmean_max =  1
		levels = np.arange(-1,1.1,0.2)

	if   varname not in ['Prec1', 'CC1']: cmap = plt.get_cmap('RdBu_r')
	elif varname     in ['Prec1', 'CC1']: cmap = plt.get_cmap('RdBu')

	print('levels -->', levels)

	print()
	print('# Plotting maps')
	print()
	
	dir_out = './figures/{}'.format(varname)
	
	if not os.path.exists(dir_out): os.makedirs(dir_out)
	
	# clevs = np.linspace(var_timmean_min, var_timmean_max, 20)
	clevs = levels
		
	projection = ccrs.Robinson(central_longitude=0)
	
	transform  = ccrs.PlateCarree(central_longitude=0)

	for time_step in range(0, len(tim), freq):
		 
		# varp = var[time_step,:,:]

		figname  = '{}/aeolus2_exps_{}_{:03}.png'.format(dir_out, varname, time_step)
		
		figtitle = '{} ({}) \ntime step {:03}'.format(varname.upper(), longname, time_step)

		# if os.path.exists(figname): continue

		# Set the figure size, projection, and extent
		
		fig = plt.figure(figsize=(8,4))

		ax = plt.axes(projection=projection)
		ax.set_global()
		ax.coastlines(resolution="110m", linewidth=1)
		ax.gridlines(linestyle='--', color='gray') #, draw_labels=True, xlocs=np.arange(-150, 151, 60), ylocs=np.arange(-90, 91, 15))

		# Set contour levels, then draw the plot and a colorbar
		
		# clevs = np.linspace(var_timmean_min, var_timmean_max, 20)

		plt.contourf(lon, lat, var[time_step,:,:], clevs, transform=transform, cmap=cmap, extend='both')

		plt.title(figtitle, size=14)
		cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8, format='%.3f')
		cb.set_label('[{}]'.format(unit), size=10, rotation=-90, labelpad=15)
		cb.ax.tick_params(labelsize=10)
		
		plt.tight_layout()
		fig.savefig(figname, format='png', dpi=60, bbox_inches='tight')
		plt.close()

		if os.path.exists(figname): print('done -->', figname)

	# exit()

	# Animations
	
	fignames = '{}/aeolus2_exps_{}_*.png'.format(dir_out, varname)
	gifname  = '{}/aeolus2_exps_{}.gif'  .format(dir_out, varname)
	mp4name  = '{}/aeolus2_exps_{}.mp4'  .format(dir_out, varname)
	
	if False:
			
		print()
		print('# Creating GIF animation from PNGs')
		print()
		
		# Create the frames
		frames = []
		imgs = sorted(glob.glob(fignames))
		
		new_frame = Image.open(imgs[0])
		width, height = new_frame.size
		
		frames.append(new_frame)
		frames.append(new_frame)
		for i in imgs:
			new_frame = Image.open(i)
			frames.append(new_frame)
		frames.append(new_frame)
		frames.append(new_frame)

		# Save into a GIF file that loops forever
		frames[0].save(gifname, format='GIF', append_images=frames[1:], save_all=True, duration=600, loop=0)

		# os.system('rm Aeolus2_exps_{}_{}.png'.format(varname, '*'))
		
		if os.path.exists(gifname): print('done -->', gifname)
		
	if False:
		
		print()
		print('# Creating MP4 animation from GIF', int(width*0.9), int(height*0.9))
		print()
				
		clip = mp.VideoFileClip(gifname, audio=False, target_resolution=(int(height*0.9), int(width*0.9)))
		clip.write_videofile(mp4name)
		
		if os.path.exists(mp4name): print('done -->', mp4name)
		
dat.close()


def dummy_image_with_text(w, h, text='text', fontsize=24, text_color='black', background='white'):
	
	blank_image = Image.new('RGBA', (w, h), background)
	img_draw = ImageDraw.Draw(blank_image)
	font = ImageFont.truetype("DejaVuSerif.ttf", fontsize)
	# img_draw.rectangle((70, 50, 270, 200), outline='red', fill='blue')
		
	textwidth, textheight = int(len(text)*fontsize*0.4), fontsize
	margin = int(fontsize/2)
	x = int(w/2) - int(textwidth/2) - int(margin/2)
	y = int(h/2) - margin #- int(textheight/2)
	
	img_draw.text((x, y), text, fill=text_color, font=font)
	# blank_image.save('drawn_image.png')
	return blank_image

# dummy_image = dummy_image_with_text(500, 50, text='text here', fontsize=20)


def pil_grid(images, max_horiz=np.iinfo(int).max, title='Title here for the figures         '):
	
	images = [Image.open(x) for x in images]
	n_images = len(images)
	n_horiz = min(n_images, max_horiz)
	h_sizes, v_sizes = [0] * n_horiz, [0] * (n_images // n_horiz)

	image_widch = images[0].size[0]

	dummy_image = dummy_image_with_text(image_widch*max_horiz, 70, text=title, fontsize=30)

	for i, im in enumerate(images):
		h, v = i % n_horiz, i // n_horiz
		h_sizes[h] = max(h_sizes[h], im.size[0])
		v_sizes[v] = max(v_sizes[v], im.size[1])
	h_sizes, v_sizes = np.cumsum([0] + h_sizes), np.cumsum([0] + v_sizes)
	
	im_grid = Image.new('RGB', (h_sizes[-1], v_sizes[-1]+dummy_image.size[1]), color='white')
	
	im_grid.paste(dummy_image, (0, 0))
	
	for i, im in enumerate(images):
				
		im_grid.paste(im, (h_sizes[i % n_horiz], v_sizes[i // n_horiz]+dummy_image.size[1]))
		
	return im_grid




# ['h1', 'u1ph', 'u1th', 'h2', 'RT1', 'RT2', 'u2ph', 'u2th', 'q1', 'b1', 'b2', 'q2', 'w1', 'w2', 'menthalpy', 'CC1', 'Prec1', 'DD1', 'DDr', 'Ev']

seq_1 = [  'b1',   'q1',    'h1',  
		   'b2',   'q2',    'h2' ]
seq_2 = ['u1ph', 'u1th', 'Prec1', 
		 'u2ph', 'u2th',   'CC1' ]
seq = seq_1 + seq_2


print()
print()
print('# Variables loop - painel with all variables in one step')
print()

for time_step in range(0, len(tim), freq):
	
	images = []
	
	for varname in seq:
		
		dir_out = './figures/{}'.format(varname)
		figname  = '{}/aeolus2_exps_{}_{:03}.png'.format(dir_out, varname, time_step)
		
		images.append(figname)
		
	# print(images)
	
	dir_out_painel = './figures/painel'
	
	if not os.path.exists(dir_out_painel): os.makedirs(dir_out_painel)
	
	painel_name  = '{}/aeolus2_exps_variables_painel_{:03}.png'.format(dir_out_painel, time_step)
	
	# if os.path.exists(painel_name): continue
	
	figtitle = 'Aeolus2.0 Experiment - {}{}'.format(sufix.upper().replace('_',' '), ' '*45)
	
	images_combined = pil_grid(images, max_horiz=3, title=figtitle)
	images_combined.save(painel_name)
	images_combined.close()
	
	if os.path.exists(painel_name): print('done -->', painel_name)



# Animations


gifname       = '{}/aeolus2_exps_variables_painel.gif'      .format(dir_out_painel)
gifname_small = '{}/aeolus2_exps_variables_painel_small.gif'.format(dir_out_painel)
mp4name       = '{}/aeolus2_exps_variables_painel.mp4'      .format(dir_out_painel)
movname       = '{}/aeolus2_exps_variables_painel.mov'      .format(dir_out_painel)


painel_fnames = []
for time_step in range(0, len(tim), freq):			
	painel_fname  = '{}/aeolus2_exps_variables_painel_{:03}.png'.format(dir_out_painel, time_step)
	painel_fnames.append(painel_fname)
		

if 1:
		
	print()
	print('# Creating GIF animation from PNGs')
	print()
	
	# Create the frames
	frames = []
	imgs = painel_fnames
	
	new_frame = Image.open(imgs[0])
	width, height = new_frame.size
	
	frames.append(new_frame)
	frames.append(new_frame)
	for i in imgs:
		new_frame = Image.open(i)
		frames.append(new_frame)
	frames.append(new_frame)
	frames.append(new_frame)

	# Save into a GIF file that loops forever
	frames[0].save(gifname, format='GIF', append_images=frames[1:], save_all=True, duration=400, loop=0)

	# os.system('rm Aeolus2_exps_{}_{}.png'.format(varname, '*'))
	
	if os.path.exists(gifname): print('done -->', gifname)

if 1:
	
	print()
	print('# Creating GIF SMALL animation from PNGs')
	print()


	# Create the frames
	frames = []
	imgs = painel_fnames
	
	new_frame = Image.open(imgs[0])
	width, height = new_frame.size
	new_width  = int( width*0.4)
	new_height = int(height*0.4)
	
	print('resolution -->', int(new_width), int(new_height))
	
	for i in imgs:
		new_frame = Image.open(i)
		new_frame = new_frame.resize((new_width, new_height))
		frames.append(new_frame)
	frames.append(new_frame)
	frames.append(new_frame)

	# Save into a GIF file that loops forever
	frames[0].save(gifname_small, format='GIF', append_images=frames[1:], save_all=True, duration=400, loop=100, optimize=True)
	
	new_frame.close()
	frames = None
	
	if os.path.exists(gifname_small): print('done -->', gifname_small)


if 0:
	
	print()
	print('# Creating MP4 animation from GIF', int(width*0.9), int(height*0.9))
	print()
			
	clip = mp.VideoFileClip(gifname, audio=False, target_resolution=(int(height*0.5), int(width*0.5)))
	clip.write_videofile(mp4name)
	
	if os.path.exists(mp4name): print('done -->', mp4name)
	

if 0:
	
	print()
	print('# Creating MOV animation from GIF', int(width*0.5), int(height*0.5))
	print()
			
	clip = mp.VideoFileClip(gifname, audio=False, target_resolution=(int(height*0.5), int(width*0.5)))
	clip.write_videofile(movname, codec='libx264')	 # this also works for mp4, but mov plays in whatsapp 
	clip.close()
	
	if os.path.exists(movname): print('done -->', movname)
	
	
exit()


















