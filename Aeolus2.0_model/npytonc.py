"""
npytonc.py
Convert a series of .npz files into one NetCDF file

Stefan Petri petri@pik-potsdam.de
"""

import sys
import os
import numpy as np
from NetCDFOutput import NetCDFOutput
from glob import glob

def usage(myname):
    print('usage: ', myname, '<NCoutfile> npz-infiles ...')
    quit()

# if (len(sys.argv) <= 2): usage(sys.argv[0])
# ncFile = sys.argv[1]
# npzFiles = sys.argv[2:]

ncFile   = 'output/teste_npytonc_output.nc'
npzFiles = glob('output/output_1*.npz')

print(npzFiles)
# exit()

if os.path.exists(ncFile):
    print("error: ",ncFile," already exists, will not overwrite it")
    usage(sys.argv[0])

if not os.path.exists(npzFiles[0]):
    print("error: cannot read input file ",npzFiles[0])
    usage(sys.argv[0])

# open first input file to obtain axis information
d = np.load(npzFiles[0])
lons = d['lamda']/(2.*np.pi)*360.0 - 180.0
phi = np.pi/2.0-d['theta']
lats = phi/(2.*np.pi)*360.0

# now create the NetCDF output file
gridinfo = {
   "lons"         : lons,
   "lats"         : lats,
}

print('create ', ncFile)
ncout = NetCDFOutput(comm=None, coords=gridinfo, f=ncFile, isrestart=False,
                     comment="from npytonc.py")
# create metadata
for v in d.files:
    if (v == 't'): continue
    if (v == 'lamda'): continue
    if (v == 'theta'): continue
    print('AddVariable ',v)
    if (d[v].ndim == 0):
        ncout.AddVariable(v, axes=(ncout.timeaxisname))
    elif (d[v].ndim == 2):
        ncout.AddVariable(v, axes=(ncout.timeaxisname,ncout.lataxisname,ncout.lonaxisname,))
    else:
        print('Warning: cannot handle variable ',v,' with ', d[v].ndim,' dimensions')

# now copy the data itself
for f in npzFiles:
    print('load ', f)
    d = np.load(f)
    for v in d.files:
        if (v == 't'): continue
        if (v == 'lamda'): continue
        if (v == 'theta'): continue
        ncout.WriteVar(v, d[v], d['t'])

ncout.close
print('generated ', ncFile)
