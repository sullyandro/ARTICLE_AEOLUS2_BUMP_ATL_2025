"""
NetCDFOutput.py

Stefan Petri petri@pik-potsdam.de
"""

#import os
import time
import numpy as np
#import dedalus.public as de
from mpi4py import MPI
from netCDF4 import Dataset
#import cftime


def make_output_filename(filename, comm):
    """
Check given filename. Append .nc if not already present.
In parallel runs, append MPI rank number in the given MPI communicator, zero-filled to 4 digits
    """
    if (filename[-3:] != ".nc"):
        filename += ".nc"
    if (comm):
        if (comm.size > 1):
            #filename += ".%(rank)04d" % { "rank" : comm.size - comm.rank - 1}
            filename += ".%(rank)04d" % { "rank" : comm.rank}
    return filename
# end of make_output_filename()

class NetCDFOutput:
    """
An object of class NetCDFOutput represents one NetCDF output file,
together with all the Variables that are to be written into it.
First, Variables have to be registered for a NetCDFOutput instance.
At that time, the corresponding metadata is created.
Afterwards, the content values of all registered Variable are written
with each invocation of OutputToFile().

For diagnostics, float is sufficient and saves half the space compared
to double precision, but for restart files, we want unconverted
doubles, in the assumption that computing across restarts should give
bitwise identical results compared to simulation runs without restart
inbetween.  However, for axes coordinates, single precision should
always be sufficient, because those are not used read from restart
files.

For MPI-parallel program runs, each MPI task writes the data of its
local domain to a separate NetCDF file, which gets the task number as
name suffix. The per-domain files can then be combined into one global
file using the mppnccombine postprocessing tool from FMS.

There is a test program to exercise and debug class NetCDFOutput in
sequential and in MPI-parallel mode at the bottom of
NetCDFOutput.py. Also, there are instructions how to compile and run
the test program.
    """

    def __init__(self, comm, coords, f,
                 isrestart=False,
                 comment="",
                 timeunits="days since 0001-01-01 00:00",
                 calendar=None):
        """
param comm is an MPI communicator, or None for sequential runs.
       There is no default value for this parameter,
        because we dont want the user to forget passing it.
param coords is a dictionary that can contain grid information:
        lons,lats,levels,lonb,latb,levelb,ics,ice,jcs,jce
param f is the filename.
 A suffix ".nc" is appended to the filename f, if it is not already there.
param isrestart specifies wether this is a restart file (true) or a diagnostic output (false).
 Restart files are written in double precision, while single floats are used for diagnostics.
param timeunits defaults the relative time axis
param comment is written to the NetCDF file as global attribute 'creator'.
 It should contain some information about program version and kind of experiment.

#vars is a dict of variables in this NetCDF file. ['varname':ncVar]
variables: The variables dictionary maps the names of variables defined for this
Dataset to instances of the Variable class.
The attributes name, long_name, units, missing_value, _FillValue, and axes
are assigned to the NetCDF variable object denoted by ncVar in AddVariable().
        """
        #self.vars = dict()
        self.record_counter = 0 # number of time steps in the NC file
        self.comm = comm
        self.isrestartfile = isrestart
        self.vardtype = 'float64' if self.isrestartfile else 'float32'
        self.filename = make_output_filename(f, comm)
        # Sigh. mppnccombine segfaults on files with default format="NETCDF4"
        # However, it works fine with format="NETCDF4_CLASSIC"
        self.datafile = Dataset(self.filename, 'w', parallel=False, format="NETCDF4_CLASSIC")
        self.datafile.history = time.ctime(time.time()) + ": Created by Aeolus2"
        self.datafile.Conventions = "CF-1.0"
        self.lonaxisname = "longitude"
        self.lataxisname = "latitude"
        self.levaxisname = "level"
        self.timeaxisname = "time"
        if (comment):
            self.datafile.creator = comment
        if (comm and comm.size > 1):
            # Each MPI task writes data only for its local compute domain.
            # Halos are never written.
            # We use the postprocessing tool mppnccombine from FMS/MOM4 to combine the per-task
            # files into one global file.
            self.datafile.NumFilesInSet = comm.size

        if (0 == coords["lons"].size):
            print("Error: axis ",self.lonaxisname,"has length 0")
        lonDim = self.datafile.createDimension(self.lonaxisname, coords["lons"].size)
        if ("lonb" in coords):
            if (0 == coords["lonb"].size):
                print("Error: axis ","lonb","has length 0")
            lonbDim = self.datafile.createDimension("lonb", coords["lonb"].size)
        if (0 == coords["lats"].size):
            print("Error: axis ",self.lataxisname,"has length 0")
        self.latDim = self.datafile.createDimension(self.lataxisname, coords["lats"].size)
        if ("latb" in coords):
            if (0 == coords["latb"].size):
                print("Error: axis ","latb","has length 0")
            latbDim = self.datafile.createDimension("latb", coords["latb"].size)
        if ("levels" in coords):
            if (0 == coords["levels"].size):
                print("Error: axis ",self.levaxisname,"has length 0")
            levelsDim = self.datafile.createDimension(self.levaxisname, coords["levels"].size)
        if ("levelb" in coords):
            if (0 == coords["lons"].size):
                print("Error: axis ","levelb","has length 0")
            levelbDim = self.datafile.createDimension("levelb", coords["levelb"].size)
        # record dimension without fixed size
        timeDim = self.datafile.createDimension(self.timeaxisname, None)

        # Define the coordinate variables
        # Sigh. CF defines attribute axis to identify the directions,
        # but the MOM tools use attribute cartesian_axis
        # Note that we _always_ use ncFloat for the space axes
        lonVar = self.datafile.createVariable(self.lonaxisname, 'f4', (self.lonaxisname,))
        lonVar.units = "degrees_east"
        lonVar.axis = "X"
        lonVar.cartesian_axis = "X"
        if ("lonb" in coords):
            lonVar.edges = "lonb"
            lonbVar = self.datafile.createVariable("lonb", 'f4', ("lonb",))
            lonbVar.units = "degrees_east"
        latVar = self.datafile.createVariable(self.lataxisname, 'f4', (self.lataxisname,))
        latVar.units = "degrees_north"
        latVar.axis = "Y"
        latVar.cartesian_axis = "Y"
        if (comm and comm.size > 1):
            # "domain_decomposition = #0, #1, #2, #3" attribute
            # #0 starting position of original dimension
            # #1 ending position of original dimension
            # #2 starting position of decomposed dimension
            # #3 ending position of decomposed dimension
            latVar.domain_decomposition = (1, coords["nlats_global"], coords["jcs"]+1, coords["jce"])
        if ("latb" in coords):
            latVar.edges = "latb"
            latbVar = self.datafile.createVariable("latb", 'f4', ("latb",))
            latbVar.units = "degrees_north"
            if (comm and comm.size > 1):
                latbVar.domain_decomposition = (1, coords["nlats_global"]+1, coords["jcs"]+1, coords["jce"]+1)
        if ("levels" in coords):
            levelsVar = self.datafile.createVariable(self.levaxisname, 'f4', (self.levaxisname,))
            levelsVar.units = "m"
            levelsVar.axis = "Z"
            levelsVar.cartesian_axis = "Z"
        if ("levelb" in coords):
            levelsVar.edges = "levelb"
            levelbVar = self.datafile.createVariable("levelb", 'f4', ("levelb",))
            levelbVar.units = "m"
        # precision of float is insufficient to store seconds as fraction of day
        self.timeVar = self.datafile.createVariable(self.timeaxisname, 'f8', (self.timeaxisname,))
        self.timeVar.units = timeunits
        self.timeVar.axis = "T"
        self.timeVar.cartesian_axis = "T"
        if (calendar):
            self.timeVar.calendar = calendar
            self.timeVar.calendar_type = calendar

        # now write out the axis data
        #print('coords["lons"]', coords["lons"])
        lonVar[:] = coords["lons"]
        if ("lonb" in coords):
            lonbVar[:] = coords["lonb"]
        #print('coords["lats"]', coords["lats"])
        latVar[:] = coords["lats"]
        if ("latb" in coords):
            latbVar[:] = coords["latb"]
        if ("levels" in coords):
            levelsVar[:] = coords["levels"]
        if ("levelb" in coords):
            levelbVar[:] = coords["levelb"]

        self.datafile.sync()
    # end of __init__()


    def AddVariable(self, var, longname="", units="", miss_val=None, axes=None):
        """Add one variable to a NetCDFOutput instance.

        var is the variable name as string (TODO: handle Dedalus Fields ?)
        axes is a tuple of axis names ("time", "level", "latitude", "longitude",)
        """
        if (var in self.datafile.variables):
            print("Warning: Variable ", var, " was already added to NetCDFOutput ", self.filename)
            print(" -- ignoring repeated AddVariable")

        ncvar = self.datafile.createVariable(var, self.vardtype,
                                     dimensions=axes,
                                     # zlib, complevel, shuffle,
                                     fill_value = miss_val)
        if (longname):
            ncvar.long_name = longname
        if (units):
            ncvar.units = units
        if (miss_val):
            # not allowed here, already done inside createVariable()
            #ncvar._FillValue = miss_val
            ncvar.missing_value = miss_val

        #self.vars[var] = ncvar
        self.datafile.sync()
    # end of AddVariable()


    def WriteVar(self, varname, data, timestamp=None):
        """timestamp can be None for static variables"""
        if (not (varname in self.datafile.variables)):
            print("Error: Variable ", varname, " was not added to NetCDFOutput ", self.filename)
            raise Exception

        ncvar = self.datafile.variables[varname]

        ## Sigh. transpose() creates output upside down.
        #print("WriteVar(",varname,",[..],",timestamp,")  dtype", data.dtype, " shape ", data.shape,
        #      " to type ",ncvar.datatype, " shape ",ncvar.shape)
        if ("time" in ncvar.dimensions):
            # find timeslice
            if (self.record_counter > 0 and timestamp < self.timeVar[self.record_counter-1]):
                print("Error: cannot write ", varname, " before last time slice: ",
                      timestamp, " < ", self.timeVar[self.record_counter-1])
                raise Exception
            if (self.record_counter == 0 or timestamp > self.timeVar[self.record_counter-1]):
                self.timeVar[self.record_counter] = timestamp
                self.record_counter += 1

            # write data
            # TODO: check if we can abbreviate the switching on ndims by using
            #       elipsisis slices [...] ?
            #print("  at time record ",self.record_counter-1)
            #if ('u1th' == varname):
            #    print("    min ",np.min(data)," max ",np.max(data))
            if (ncvar.ndim == 0): # static scalar
                print("cannot write static scalar to time-dependent variable ", varname)
                raise Exception
            elif (ncvar.ndim == 1): # timeseries scalar (presumably)
                if isinstance(data, np.ndarray):
                    ncvar[self.record_counter-1] = data.astype(self.vardtype)
                else:
                    ncvar[self.record_counter-1] = data
            elif (ncvar.ndim == 2): # MeridionalNcvariable1D T,Y
                ncvar[self.record_counter-1,:] = data[:].astype(self.vardtype)
            elif (ncvar.ndim == 3): # HorizontalNcvariale T,Y,X or MeridionalNcvariable T,Z,Y
                # TODO: compare shapes to find need/possibility for transposition
                foo = np.transpose(data)[:   ,:].astype(self.vardtype)
                #print("  foo dtype", foo.dtype, " shape ", foo.shape)
                ncvar[self.record_counter-1,:,:] = foo
            elif (ncvar.ndim == 4): # 3-D Ncvariable T,Z,Y,X
                ncvar[self.record_counter-1,:,:,:] = np.transpose(data)[:,:   ,:].astype(self.vardtype)
            else:
                raise Exception # too many dimensions

        else: # write static data
            #print("  as time-independent static data ")
            if (ncvar.ndim == 0): # static scalar
                #try:
                #    ncvar[...] = data.astype(self.vardtype)
                #except AttributeError:
                #    ncvar[...] = data
                if isinstance(data, np.ndarray):
                    ncvar[...] = data.astype(self.vardtype)
                else:
                    ncvar[...] = data
            elif (ncvar.ndim == 1): # MeridionalNcvariable1D Y
                ncvar[:] = data[:].astype(self.vardtype)
            elif (ncvar.ndim == 2): # HorizontalNcvariale Y,X or MeridionalNcvariable Z,Y
                #print("data.shape ", data.shape)
                # TODO: compare shapes to find need/possibility for transposition
                ncvar[:,:] = np.transpose(data)[:   ,:].astype(self.vardtype)
            elif (ncvar.ndim == 3): # 3-D Variable Z,Y,X
                ncvar[:,:,:] = np.transpose(data)[:,:   ,:].astype(self.vardtype)
            else:
                raise Exception # too many dimensions

        self.datafile.sync()

    # end of WriteVar()

    def close(self):
        """
        Just as the function name says.
        """
        self.datafile.close()
    # end of close()

# end of class NetCDFOutput


def test_NetCDFOutput():
    """
    (unset I_MPI_DAPL_UD I_MPI_PMI_LIBRARY I_MPI_DAPL_UD_PROVIDER; python3 -u Aeolus2_1.py -i 10 -o 5) >& o

    (unset I_MPI_DAPL_UD I_MPI_PMI_LIBRARY I_MPI_DAPL_UD_PROVIDER; python3 -u -c 'from NetCDFOutput import NetCDFOutput, test_NetCDFOutput; test_NetCDFOutput()')

    """

    print("testing NetCDFOutput()")
    nlons = 96
    nlats = 48

    comm      = MPI.COMM_WORLD
    rank      = comm.rank
    size      = comm.size
    print("MPI rank ", rank, " of ", size)

    #dlon = 360/nlons
    #dlat = 180/nlats
    ics = 0
    ice = nlons
    lons = np.linspace(0,360,nlons,dtype=float)
    #jcs = 0 if rank > 0 else int(nlats/rank)
    if rank > 0:
        jcs = int(nlats/rank)
    else:
        jcs = 0
    jce = int(nlats/(rank+1))
    lats = np.linspace(-90,90,nlats,dtype=float)[jcs:jce]

    gridinfo = {
        "nlons_global" : nlons,
        "nlats_global" : nlats,
        "lons" : lons,
        "lats" : lats,
        "ics"   : ics,
        "ice"   : ice, # according to python conventions, one after end
        "jcs"   : jcs,
        "jce"   : jce # according to python conventions, one after end
    }
    print(gridinfo)

    ncout = NetCDFOutput(comm, gridinfo, "test_Aeolus2-output.nc",
                     comment="Testing NetCDFOutput.py")

    # Here we must use the names of the axes as defined in the NetCDF file
    # TODO? implement a function to retrieve the axis names.
    axes_tyx = ("time","latitude","longitude",)
    axes_yx = ("latitude","longitude",)

    ncout.AddVariable('test1', 'long name of test1', 'units of test1', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('test2', 'long name of test2', 'units of test2', miss_val=None, axes=axes_yx)

    t = 42.42

    ncout.WriteVar('test1', np.array([[1]]), 1.0)
    ncout.WriteVar('test1', np.array([[2]]), 2.0)
    ncout.close()
    return

    #test1 = np.indices((ice-ics, jce-jcs))
    #print(test1)
    #ncout.WriteVar('test1', test1, t)
    #ncout.WriteVar('test2', test1, timestamp=None)

# end test_NetCDFOutput()
