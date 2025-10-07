"""
Aeolus2

The dynamical core of this atmosphere model is based on pseudo-spectral
moist-convective Thermal Rotating Shallow Water (mcTRSW) model on the
full sphere. Numerical methods of the model is according to Dedalus
algorithm that employs spin-weighted spherical harmonics.
Detailed explanation of the model is in the following article:
Rostami, M., Zhao, B. and Petri, S. (2022), On the genesis and dynamics
of MJO-like structure formed by equatorial adjustment of localized heating.
Q J R Meteorol Soc. 148(749) pg 3788--3813. https://doi.org/10.1002/qj.4388

The functions provided in this code are intended to be invoked from a coupler
interface, as well as from a stand-alone main program.
"""

import sys
import os
import time
start_time = time.time()
#import argparse
#import pathlib
import math                     as ma
#import mpmath                   as mp

import numpy                    as np
import scipy.integrate          as integrate #Gives access to the ODE integration package
import scipy.special            as ss
import scipy.io                 as sio
#from scipy.sparse import linalg as spla

from mpi4py import MPI
#from multiprocessing import Pool

from netCDF4 import Dataset # for restart files

import dedalus.public           as de
# Load config options
#from dedalus.tools.config import config

import sphere_wrapper           as sph
import equations_aeolus2main    as eq
import timesteppers
from NetCDFOutput import NetCDFOutput

# Import user settings from config.py file
# TODO: avoid the wildcard import, for better readabilty and maintainability
#print('sys.path', sys.path)
# Try to increase the probability that config.py is actually taken from
# the current experiment directory and not from the source code directory.
sys.path.insert(0, os.getcwd())
from config import *            # import constants and parameters from config file
import  config as config_ael
varconf = vars(config_ael) # do this on global scope

#--------- *** insolation and the orbital parameters over the last 5 million years *** -----------
import climlab
from   climlab                  import constants as const
from   climlab.solar.insolation import daily_insolation
from   climlab.solar.orbital    import OrbitalTable
from   climlab.utils.thermo     import pseudoadiabat
#----------------------------------------

epsilon = 1e-10
def double_eq(v1,v2,e=epsilon): return np.fabs((v1-v2)/((v1+v2)/2)) < e

# In both stand-alone and coupled configurations,
# this module does communication only in its own communicator,
# which is not necessarily COMM_WORLD.
# If coupled to FMS, the invoking aeolus2fms module takes care
# of setting up the inter-communicator for exchange with FMS.
comm      = None
rank      = None
size      = None

# Define some physical constants
R_dry             = 287.05   # R_dry is the specific gas constant for dry air [J/(kg*K)].
R_moist           = 461.50   # R_moist​ is the specific gas constant for water vapor [J/(kg*K)].
tau_bu            = 10.0     # numerical relaxation time for smooth buoyancy forcing



## Define some variables on global scope.
L_dealias     = None
L_max         = None # spherical harmonic order
S_max         = None # spin order (leave fixed)
domain        = None
S             = None # Sphere object
Nxx           = None # number of longitudes in global grid
Nyy           = None # number of latitudes in global grid
ni            = None # number of longitudes in local grid
nj            = None # number of latitudes in local grid

# indices of the borders of the compute domain in the global grid.
# Note: following python conventions, the end indices are one
# _beyond_ the actual end.
ics  = 0 # sigh. "is" is a reserved word in python
ice  = None
jcs  = 0
jce  = None
# longitudes and latitudes for the local compute domain
lons_local = None
lats_local = None


ncout = None     # NetCDFOutput handle
axes_t = None    # time series of scalars
axes_ty = None   # time series of zonal-mean values
axes_yx = None   # static maps
axes_tyx = None  # time series of maps

HORO = None      # height of orography [m]
SIGORO = None    # sub-grid standard deviation of orography height
land_frac = None # fraction of land in grid cell, Ocean = 1-land_frac
# cell area and area_sum are needed for properly weighting global mean values.
cell_area = 1.0  # for e.g. properly weighted global mean values
area_sum = 1.0   # global sum of cell areas, i.e. planet surface area

# start time of experiment
# FMS time stamps are broken out into days and seconds
Time_init_year = None
Time_init_month = None
Time_init_day = None
Time_init_seconds = None
# start time of this run, within the experiment
#Time_days = None
#Time_seconds = None

Time_step_days = None
Time_step_seconds = None

without_topo = False
update_land_frac_always = False

nr_tracers = 0  # number of tracers to exchange with coupler
q_ind = -1      # index of humidity tracer
co2_ind = -1    # index of CO2 tracer

# Tensor variables which are used with the Dedalus engine to
# compute / advance the state of the model.
u1   = None
h1   = None
u2   = None
h2   = None
q1   = None
b1   = None
b2   = None
q2   = None
w1   = None
w2   = None

# some more Tensor variables.
# TODO: check which of these could be replaced by 2-D numpy arrays
# instead of the more expensive Tensor structures.
u10  = None # remains constant after init
u20  = None # remains constant after init
b10  = None # remains constant after init
b20  = None # remains constant after init
h10  = None # accessed in coefficient space in nonlinear()

Du1  = None # accessed in coefficient space in nonlinear()
Du2  = None # accessed in coefficient space in nonlinear()
Db1  = None # accessed in coefficient space in nonlinear()
Db2  = None # accessed in coefficient space in nonlinear()
Dh1  = None # accessed in coefficient space in nonlinear()
Dh10 = None # accessed in coefficient space in nonlinear()
Dh2  = None # accessed in coefficient space in nonlinear()

uh1  = None
uh2  = None
uq1  = None
uq2  = None
#ub1  = None
uhf1 = None
uqf1 = None
uqf2 = None
ubf1 = None
uw1  = None
uw2  = None
uwf1 = None
uwf2 = None

divuh1  = None
divuh2  = None
divuq1  = None
divuq2  = None
divuw1  = None
divuw2  = None
#divub1  = None
#divuhf1 = None
#divuqf1 = None
#divuwf2 = None
#divubf1 = None


u1_rhs = None
h1_rhs = None
u2_rhs = None
h2_rhs = None
q1_rhs = None
b1_rhs = None
b2_rhs = None
q2_rhs = None
w1_rhs = None
w2_rhs = None

Ev     = None # evaporation
Oc     = None # Ocean
CC1    = None # Condensation/Latent heat release forcing/ep
DD1    = None # Downdraft+bias correction due to condensation in the lower layer, from upper troposphere to the lower troposphere
DDr    = None # Downdraft without bias correction due to condensation in the lower layer, from upper troposphere to the lower troposphere
VV2    = None # Vaporization in the upper layer
VV1    = None # Vaporization in the lower layer
Prec1  = None # Precipitation in the lower layer
speed  = None # velocity magnitude
R2     = None
Rad1   = None # Newtonian cooling for the lower layer
Rad2   = None # Newtonian cooling for the upper layer
RT1    = None # RRTM_LW/SW, lower layer
RT2    = None # RRTM_LW/SW, upper layer
q_o_h1 = None # q1/h1
q_o_h2 = None

# arrays for pressure at layers 1 and 2
Press1 = None
Press2 = None

state_vector = None
RHS          = None
timestepper  = None

step = None
dt = None
t = 0 # timestep

# build matrices
P,M,L,LU = [],[],[],[]


def aeolus2_init_grid(nlons_global, nlats_global, lon_0_360, comm_in):
    """
    Initialise the Dedalus grid accroding to the incoming
    parameters.
    Set global variables with information about communicator,
    domain, boundary indices, etc.
    Also initializes the NetCDF output file.
    Return domain boundary indices, and axis labels in degrees.
    """

    global L_max, L_dealias, S_max, domain, S, Nxx, Nyy
    global comm, rank, size
    global ics, ice, jcs, jce
    global ni, nj
    global ncout, axes_t, axes_ty, axes_yx, axes_tyx
    global step

    L_dealias = 3/2
    S_max     = 3    # spin order (leave fixed)

    comm      = comm_in
    rank      = comm.rank
    size      = comm.size

    # print some information about the execution environment
    if 0 == rank:
        print('python version ',sys.version)
        print('python installed in ',sys.executable)
        print('running on ',sys.platform)
        print('python path is ',sys.path)
        print('float info', sys.float_info)
        print('Numpy version ',np.version.full_version)
        print('Numpy default floating point error handling ',np.geterr())
        # Sigh. Dedalus provides version info only in its main module,
        # but that must not be imported here.
        # dedalus.__version__
        # dedalus.__path__
        print('Dedalus public module info ',de)
    #np.seterr(all='raise')

    if (nlons_global == 768 and nlats_global == 384):
        print("running on smooth grid ", nlons_global, nlats_global, flush=True)
        L_max = 255  # spherical harmonic order
        step = 1.   # minutes, numerical time step
    elif (nlons_global == 384 and nlats_global == 192):
        print("running on fast grid ", nlons_global, nlats_global, flush=True)
        L_max = 127  # spherical harmonic order
        step = 2.5  # minutes, numerical time step, the model is convergent at step=5
    elif (nlons_global == 192 and nlats_global == 96):
        print("running on superfast grid ", nlons_global, nlats_global, flush=True)
        L_max = 63  # spherical harmonic order
        step = 5.0 # 5 is ok., tested, numerical time step
    else:
        print("Error: dont know how to handle grid resolution ",
              nlons_global,"x",nlats_global, flush=True)
        raise Exception

    Nxx = nlons_global
    Nyy = nlats_global

    theta_basis = de.Fourier('theta', L_max+1, interval=(0,np.pi), dealias=L_dealias)
    #print('theta_basis', theta_basis)

    # Make domain
    lamda_basis = de.Fourier('lamda', 2*(L_max+1), interval=(0,2*np.pi), dealias=L_dealias)
    domain      = de.Domain([lamda_basis,theta_basis], grid_dtype=np.float64, comm=comm)

    # set up sphere
    m_start = domain.distributor.coeff_layout.start(1)[0]
    m_len   = domain.distributor.coeff_layout.local_shape(1)[0]
    m_end   = m_start + m_len - 1
    N_theta = int((L_max+1)*L_dealias)
    print("L_max, S_max, L_dealias, m_start, m_len, m_end, N_theta",
          L_max, S_max, L_dealias, m_start, m_len, m_end, N_theta, flush=True)

    S = sph.Sphere(L_max+1,S_max,N_theta=N_theta,m_min=m_start,m_max=m_end)
    #print("S grid, weights, sin_grid, cos_grid",
    #      S, S.grid, S.weights, S.sin_grid, S.cos_grid)

    lamda        = domain.grids(L_dealias)[0]
    theta_slice  = domain.distributor.grid_layout.slices(domain.dealias)[1]
    theta_len    = domain.local_grid_shape(domain.dealias)[1]
    theta_global = S.grid
    theta        = S.grid[theta_slice].reshape((1,theta_len))
    #phi          = np.pi/2.0-theta

    #print("theta_len ",        theta_len)
    print("theta_global.size .shape ", theta_global.size, theta_global.shape)
    #print("theta_slice.start/stop/step", theta_slice.start, theta_slice.stop, theta_slice.step)
    #print("theta_slice",  theta_slice)
    print("theta_global", theta_global)
    #print("theta", theta)
    #print("theta in deg", theta/(2*np.pi)*360.0)

    # Sigh. The current input topographies are organized as -180 to 180 deg longitude,
    # but lamda runs 0 to 2pi .
    # Thus we need to rotate by -180deg to obtain a correct axis for NetCDF.
    lons = lamda[:,0]/(2.*np.pi)*360.0

    if not lon_0_360:
        print("rotating longitudes to -180..180")
        lons = lons - 180.0
    else:
        print("Warning: FMS requests longitudes 0...360")
        print("TODO: Check consistency with Aeolus2 input data", flush=True)

    # Sigh. This results in Y running N to S, while other
    # FMS-compatible atmospheres have coordinates running S to N.
    # At least this is consistent with the current input topographies.
    # ncview and ferret both revert flip the latitudes automagically,
    # with more or less warnings printed out.
    lats = theta[0,:]/(2*np.pi)*360.0 - 90.0
    print("Warning: latitudes running N to S", flush=True)
    print('lats', lats)
    lats_global = theta_global[:]/(2*np.pi)*360.0 - 90.0
    print('lats_global', lats_global)

    # Indices of the local domain boundaries in the global domain.
    # Naming convention taken from GFDLs FMS coupler.
    # i/j/k denotes lon/lat/vertical direction
    # c/d for compute resp. data domain (i.e. without/with halos)
    # s/e for start/end
    ics  = 0 # sigh. "is" is a reserved word in python
    ice  = lons.size
    jcs  = theta_slice.start
    jce  = theta_slice.stop
    # Caution: again, F90/FMS counts indices from 1, while python/Sphere counts from 0
    # The end indices, however, are the same, because F90/FMS counts to last element,
    # while python/Sphere counts 1 beyond.
    ni = ice-ics
    nj = jce-jcs
    print("ics ice jcs jce ", ics, ice, jcs, jce)
    print("number of local grid cells ni, nj: ", ni, nj)
    if not np.any(lats_global[jcs:jce] == lats[:]):
        print("Error: lats_global[jcs:jce] != lats[:]", flush=True)
        raise Exception

    #print("Lon size, Lat size: ", lons.size, lats.size)
    print("local Lons shape, Lats shape: ", lons.shape, lats.shape)
    #print("lons",    lons)
    #print("lats",    lats)

    #dlons=lons[0:lons.size-2] - lons[1:lons.size-1]
    dlon=360.0/lons.size
    print('dlon', dlon)
    lonb = np.empty(lons.size+1)
    #print("lonb.shape: ", lonb.shape)
    lonb[0] = lons[0] - dlon/2
    #print("lons slices shapes ",
    #      lons[0:lons.size-1].shape," ",
    #      lons[1:lons.size].shape)
    lonb[1:lons.size] = (lons[0:lons.size-1]+lons[1:lons.size])/2
    lonb[lonb.size-1] = lons[lons.size-1] + dlon/2
    #latb = np.empty(lats.size+1)
    #latb[0] = 90.0
    #latb[1:lats.size] = (lats[0:lats.size-1]+lats[1:lats.size])/2
    #latb[latb.size-1] = -90.0

    latb_global = np.empty(lats_global.size+1)
    latb_global[0] = 90.0
    latb_global[1:lats_global.size] = (lats_global[0:lats_global.size-1]+lats_global[1:lats_global.size])/2
    latb_global[latb_global.size-1] = -90.0
    # TODO: cut out local lonb from lonb_global - sort out correct array bounds
    latb = latb_global[theta_slice.start:theta_slice.stop+1]
    # Return longitudes and latitudes of cell boundaries and cell midpoints
    # in degrees, [and of level heights in m]
    # to be used as axes descriptions in the FMS diagnostics code.
    # The parameters are 1-D arrays, defined on the global (!) domain.
    # The returned values always describe a regular grid, regardless if the
    # atmosphere actually uses a reduced grid or whatever.

    gridinfo = {
        "nlons_global" : lons.size,
        "nlats_global" : theta_global.size,
        "lons"         : lons[:],
        "lonb"         : lonb[:],
        "lats"         : lats[:],
        "latb"         : latb[:],
        "ics"          : ics,
        "ice"          : ice,  # according to python conventions, one after end
        "jcs"          : jcs,
        "jce"          : jce,  # according to python conventions, one after end
    }

    global output_folder
    output_folder = "."
    ncout = NetCDFOutput(comm, gridinfo,
                         os.path.join(output_folder,
                                      "Aeolus2-output.nc"),
                         comment="Testing NetCDFOutput.py")

    # Here we must use the names of the axes as defined in the NetCDF file
    axes_t   = (ncout.timeaxisname,)                   # time series of scalars
    axes_ty  = (ncout.timeaxisname,ncout.lataxisname,) # time series of zonal-mean values
    axes_yx  = (ncout.lataxisname,ncout.lonaxisname,)  # static maps
    axes_tyx = (ncout.timeaxisname,ncout.lataxisname,ncout.lonaxisname,) # time series of maps
    # Sigh. The mppnccombine postprocessing tool does not work
    # for files with only static variables but a time axis.
    # Thus we must create one dummy variable with one timestep.
    # This should not be necessary once we are beyond develoment/debugging.
    ncout.AddVariable('dummy', 'dummy variable for creating a time step', '', miss_val=None, axes=axes_t)
    fortytwo = np.ndarray(42, dtype=np.float32)
    ncout.WriteVar('dummy', fortytwo, timestamp=0)

    global lons_local, lats_local
    lons_local = lons[:]
    lats_local = lats[:]
    print("shape of lons_local ",lons_local.shape," of lats_local ",lats_local.shape)

    print("aeolus2_init_grid() finished")
    # due to the stripe-shaped domain partitioning, global and local lons are identical
    return ics, ice, jcs, jce, lons, lonb, lats_global, latb_global
# end of aeolus2_init_grid()


def aeolus2_init(x):
    """
    Perform the final steps of model initialisation.
    Because FMS passes lots of parameters here,
    they are encoded into a tuple, to make python
    syntax checkers happy.
    But be careful with the variable ordering!
    """

    # FMS time stamps are broken out into days and seconds
    global Time_init_year, Time_init_month, Time_init_day, Time_init_seconds
    # current time, days since beginning of the time axis
    # (not day of year nor day of month)
    #global Time_days, Time_seconds
    global Time_step_days, Time_step_seconds
    global without_topo, update_land_frac_always
    global nr_tracers, q_ind, co2_ind
    global HORO,SIGORO
    global land_frac
    global cell_area, area_sum

    # unpack the argument tuple. Be careful with the ordering!
    Time_init_year, Time_init_month, Time_init_day, Time_init_seconds, \
    Time_days, Time_seconds, \
    Time_step_days, Time_step_seconds, \
    without_topo, update_land_frac_always, \
    nr_tracers, q_ind, co2_ind, \
    HORO, SIGORO, land_frac, cell_area \
    = x

    print('aeolus2_init:')
    print('Time_init year, month, day, seconds',
          Time_init_year, Time_init_month, Time_init_day, Time_init_seconds)
    print('Time_step days, seconds',
          Time_step_days, Time_step_seconds)
    print('without_topo', without_topo)
    print('update_land_frac_always', update_land_frac_always)
    print('nr_tracers, q_ind, co2_ind', nr_tracers, q_ind, co2_ind, flush=True)

    # Define and write out topography as static 2-D variable.
    # It might be all-zeros if without_topo is set.
    ncout.AddVariable('HORO', 'effective topography height', 'm', miss_val=None, axes=axes_yx)
    ncout.AddVariable('SIGORO', 'sub-grid standard deviation of topography height', 'm', miss_val=None, axes=axes_yx)
    ncout.AddVariable('land_frac', 'fraction of land in grid cell', '1', miss_val=None, axes=axes_yx)
    ncout.AddVariable('cell_area', 'cell area', 'm^2', miss_val=None, axes=axes_yx)
    #ncout.AddVariable('ocean', 'land-sea mask', '1', miss_val=None, axes=axes_yx)
    ncout.WriteVar('HORO',  HORO)
    ncout.WriteVar('SIGORO',  SIGORO)
    ncout.WriteVar('land_frac', land_frac)
    ncout.WriteVar('cell_area',  cell_area)
    area_sum = global_sum(comm, cell_area)
    print("global cell area sum: ", area_sum, flush=True)


    # TODO: read parameter files, init Dedalus, init submodules (LWR, SWR, ...)

    # ---------- *** User setting & parameterization *** -----------
    varconf['moist_convection']   = moist_convection
    global topography, a, T_scale, H0_dim, g_real, tau_bu
    topography = 1-without_topo
    varconf['topography']         = topography
    varconf['Radiative_Transfer'] = Radiative_Transfer
    varconf['summer_solstice']    = summer_solstice
    varconf['winter_solstice']    = winter_solstice
    varconf['unparallel_TQG_adj'] = unparallel_TQG_adj
    varconf['T_s0']               = T_s0
    varconf['CO2ppmv']            = CO2ppmv
    varconf['L_d_eq']             = L_d_eq
    varconf['B1']                 = B1
    varconf['B2']                 = B2
    #varconf['n_iterations']       = n_iterations
    varconf['delta1']             = delta1
    #print('varconf[delta1] ', delta1)

    g0                    = 1.0
    second_order_q_eff    = 1                     # 1: includes effect of q_i/h_i on the rhs of h_i and b_i equations.
    vaporization_cloud    = 1                     # 1: activates vaporization + precipitable water as clouds for the lower and upper troposphere
    Newtonian_cooling     = 1                     # 1: Newtonian cooling is active, 0: thermal relaxation of h toward initial pressure level
    uniform_DD            = 0                     # 1: uniform downdraft when downdraft is 1 (ON), 0: nonuniform downdraft
    downdraft             = 1                     # 1: active downdraft, 0: no downdraft. With real bottom topo downdraft should be 1

    if topography==1:
        a               = 1.7218         # 1.7218 corresponds to H0=10km, #a=Aspect ratio with respect to barotropic equatorial Rossby deformation radius (L_d)
    else:
        a               = 1.7            # Earth radius. Aspect ratio with respect to barotropic equatorial Rossby deformation radius, L_d=sqrt(U_bt/(beta))
                                    # U_bt=sqrt(g*H0); beta= 2*omega*cos((phi/180)*pi)/R_e; R_e: Earth's radius
    T_scale            = a/2.           # Time scale to convert to [day]:e.g.:  2*pi/(beta*L_d_eq)
    Om                 = a*1.0/2.       # Omega. It affects parametrization of Coriolis; f=2*Omega*Cos(theta) and scaling
    H0_dim             = 28000.0        # Needs to be modified, scale of (pseudo)-height (m); or, e.g. (gamma/(gamma-1))*H_s, or C_p*theta_s/g, where C_p = 1000 J/kg. K
    #CO2ppmv            = 280.
    #T_s0               = 250.0277       # 262.1009#280.
    #L_d_eq             = 4262500.       # Equatorial Rossby deformation radius= sqrt(c/beta), where c=sqrt(BgH),
    beta_real          = 2.2793**(-11.) # 2*omega*cos((phi/180)*pi)/R_e, phi=lat, R_e: Earth's radius
    #sigma_Boltzmann    = 5.671**(-8.)   # the Stefan-Boltzmann constant
    g_real              = 9.8            # gravity
#Co_Cl              = 1.25           # 1.25, amplification factor of the land-ocean warming contrast,
                                    # bigger values represent bigger contrast; indeed, Sensible Heat Flux can be merged to this parameter too.
                                    # in improved versions the contrast should be larger for drier land regions and a function of lat.
    tau_bu            = 10.0            # numerical relaxation time for smotth buoyancy forcing

    # ------------ *** End of user setting *** -------------

    period           = 2.0*np.pi/Om
    gr               = g0*(B2-1.0)

    delta2  = 1.0-delta1
    H0      = 1.0
    H1      = delta1*H0
    H2      = delta2*H0
    cg1     = np.sqrt(g0*H0)

    # Integration parameters
    # Time step length is passed as
    # Time_step_days, Time_step_seconds to this init function
    #dt           = step*60.*period/86400.  # Equal to one minute numerical time s
    global dt
    dt = (Time_step_days*86400 + Time_step_seconds)/T_scale
    print('T_scale ', T_scale)
    print('dt ', dt, flush=True)

    if topography==1:
        t_p_inv      = 1./(30.0*(dt/step)) # 1./(30.*(dt/step))  # convective adjustment/relaxation time scale, 1/(condensation relaxation time)
                                     # bigger than 30, leads to higher condensation on the equator
        gamma_NC1    = 1.0             # 2.0 coefficient of the Newtonian Cooling force (can be different from gamma and for each layer)
        gamma_NC2    = 1.0             # 2.5
        gamma_RT1    = 1.#             # 1.0001, coefficient of the Radiative Transfer flux in the lower layer (can differ from gamma_* and vary for each layer)
        gamma_RT2    = 1.              # 1.0004, coefficient of the Radiative Transfer flux in the upper layer (can differ from gamma_* and vary for each layer)
                                      # *** Notion: Values should be lower than 1 and not bigger than 1.004 that leads to an abrupt variation of h1 and h2 and
                                      # unrealistic results. Lower values in Jan. --> bigger u1 and u2
        gamma_topo   = 1.0                 # 2.1, bigger = reducing the bottom topographic effects
        nu           = step*0.00008        # 2.e-5 Newtonian viscosity/diffusivity/thermoconductivity for momentum and buoyancy equations to preclude small-scale
                                      # convective instabilities
        kappa        = 1.1*nu              # the same term for h
        coeff_RT1    = 1.0 #0.7
    else:
        t_p_inv      = 1./(30.*(dt/step))  # 1/(condensation relaxation time)
        gamma_NC1    = 1.0#2.5                 # coefficient of the Newtonian Cooling force (can be different from gamma and for each layer)
        gamma_NC2    = 1.0#3.5
        gamma_RT1    = 1.0                 # coefficient of the Radiative Transfer flux in the lower layer (can differ from gamma_* and vary for each layer)
        gamma_RT2    = 1.0                 # coefficient of the Radiative Transfer flux in the upper layer (can differ from gamma_* and vary for each layer)
                                      # gamma_RT should not be much bigger than 1.
        gamma_topo   = 1.0                 # 1.7# bigger = reducing the bottom topographic effects
        nu           = step*0.0001         # 2.e-5 Newtonian viscosity/diffusivity/thermoconductivity for momentum and buoyancy equations to preclude small-scale
                                      #  convective instabilities
        kappa        = 0.7*nu              # the same term for h
        coeff_RT1    = 1.0

    t_v_inv      = (1.0/20.0)*t_p_inv     # 1/(vaporization relaxation time)
    t_drag_inv   = 0.                     # 1/tau_u, just for TQG adjustment

    if vaporization_cloud==1:
        ivap=1.0
    else:
        ivap=0.0

    if moist_convection==1:
        i_mc      = 1.0
        if topography==1:
            Qs1       = 0.9*H1#0.01#                # physically it can be fixed to 0.9 and H_i be included in ep_i but the results will be the same
            Qs2       = 0.9*H2 #0.01#0.04
            Q01       = 0. #Qs1*(1.0-0.015)
            Q02       = 0.#Qs2*(1.0-0.015)
            max_q     = 1.10*Qs1 #5.0*(Qs1-Q01)
            gamma     = 0.2                       # Lower = more intense condensation with higher latent heat release
            ep1       = 1.0*(1.0 - gamma)          # condensation efficiency =Lv/Cp/d_theta
            t_r_b1    = 0.*1.0/(22.0*T_scale)      # Relaxation to background SST, in literature it can be even more than 25 days for SST
            t_r_b2    = 0.*1.0/(27.0*T_scale)      # Bigger relaxation leads to strengthening the insolation, bigger ---> b1>b2
            wcritical = 0.0001                     # threshold value of precipitation.
        else:
            Qs1       = 0.9*H1# 0.04               # physically it can be fixed to 0.9 and H_i be included in ep_i but the results will be the same
            Qs2       = 0.9*H2 #0.04
            Q01       = Qs1*(1-0.01)               # Close to saturation value
            Q02       = Qs2*(1-0.01)
            t_e_inv   = 15. #3.0#4.0               # for conceptual cases. Indeed, it can be time independent, but also a function of SST. 0.7 is a stable parameter
            coeff_v   = 0.1*t_e_inv                # 0.1 efficiency/effect of the surface velocity upon sea surface evaporation
            gamma     = 0.1 #0.3
            t_r_b1    = 0.*1.0/(22.0*T_scale)      # 1./10*dt#0.#      /4*60*dt # in literature is more than 25 days for SST
            t_r_b2    = 0.*1.0/(27.0*T_scale)      # 1./10*dt#0.#     /4*60*dt
            ep1       = 1.0 - gamma                # condensation efficiency =Lv/Cp/d_theta
            wcritical = 0.005                      # threshold value of precipitation.
        upward    = 0.3                           # upward draft of humidity from the lower troposphere to the upper layer
        t_r_inv   = 1.0/(130.0*T_scale)           # 1.0/(90.0*T_scale)# Thermal radiation relaxation time/ Newtonian cooling relaxation time ~ 90 days, if it relaxes to initial condition, otherwise it can be bigger.
        t_r_u1    = 0.                            # applies just for getting TQG adjustment, otherwise is zero
        t_r_u2    = t_r_u1
    elif  unparallel_TQG_adj==1:
        i_mc      = 0.
        iradN     = 0.
        Qs1       = 0.9*H1
        Qs2       = 0.9*H2
        gamma     = 0.1 #0.3
        ep1       = 1.0*(1.0 - gamma)             # condensation efficiency =Lv/Cp/d_theta
        Q01       = Qs1*(1.0-0.01)
        Q02       = Qs2*(1.0-0.01)
        t_e_inv   = 0.0
        t_r_inv   = 0.                            # restoring h1
        #--------------- *** freezing u2,b1,b2 (just three of them should be freezed to get TQG balanced state) *** --------------
                                             # The best choice is freezing u2,b1,b2, if u2 is an active layer.
        t_r_b1    = 0.*1./(4*dt)                  # freezing b1, applies just for getting TQG adjustment, otherwise is zero
        t_r_b2    = 0.*1./(4*dt)                  # freezing b2, applies just for getting TQG adjustment, otherwise is zero
        t_r_u1    = 0.*1./(4*dt)                  # freezing u1, applies just for getting TQG adjustment, otherwise is zero
        t_r_u2    = 0.*1./(4*dt)                  # freezing u2, applies just for getting TQG adjustment, otherwise is zero
        #------------------------------------------ end of freezing ------------------------------

        upward    = 0.                            # upward draft of humidity from the lower troposphere to the upper layer
        wcritical = 0.005                         # threshold value of precipitation.
    else: #dry case
        i_mc      = 0.
        Qs1       = 0.9*H1
        Qs2       = 0.9*H2
        Q01       = Qs1*(1-0.01)
        Q02       = Qs2*(1-0.3)
        t_e_inv   = 0.0
        t_r_inv   = 0.                # restoring h1
        t_r_b1    = 0. #1./(5*dt)     # 1/4*60*dt#, restoring b1
        t_r_b2    = 0. #1./(5*dt)     # 1/4*60*dt#, restoring b2
        t_r_u1    = 0. #1./(4*dt)     # applies just for getting TQG adjustment, otherwise is zero
        t_r_u2    = t_r_u1
        ep1       = 0.
        upward    = 0.               # upward draft of humidity from the lower troposphere to the upper layer
        wcritical = 0.005             # threshold value of precipitation.

    if topography==1:
        itopo  = 1.0
    else:
        itopo  = 0.0

    if second_order_q_eff==1:
        iq_o_h = 0.01#20.744 # #1.0/200.#  re-scaling for arbitrary value of: q/h[dimensional]=[L_v.g/(C_p.theta_s)].q*/h*, where * indicates nondimensional value. q_real/q_model=45
    else:
        iq_o_h = 0.0

    if Newtonian_cooling==1:
        iradN  = 1.0
        iradh  = 0.0
    else:
        iradN  = 0.0
        iradh  = 0.0 # 1.0: for choosing another cooling option
    if Radiative_Transfer==1:
        iRT    = 1.0
    else:
        iRT    = 0.0

    if downdraft==1:
        iDD    = 1.0
    else:
        iDD    = 0.0

    Qe1   = Qs1                      # saturation Q for evaporation
    Qe2   = Qs2

    #------------- *** Begining of insolation *** --------
    # The FMS coupler already provides cosz (cosine of the zenith angle)
    # and TOA solar irradiation.
    # Thus we dont need to compute that here.
    # For stand-alone Aeolus2, the daily_insolation must be computed in
    # that main program.
    # ----- *** insolation for millenial scale, present, or past *** ------ :
    #orb_0      = OrbitalTable.interp(kyear=0)                 # present-day orbital parameters
    #orb_10     = OrbitalTable.interp(kyear=-10)               # orbital parameters for 10 kyrs before present
    #days       = np.linspace(0, const.days_per_year, 365)
    #insol      = daily_insolation(lats[0,:], days, orb_0)     # lat should be from -90 to 90 degrees
    ##insol      = daily_insolation(-lats[0,:], days, orb_0)
    #if Radiative_Transfer==1:
    #    pass
    #else:
    #    insol      = 0.*insol

    #print('The area-weighted global, annual average of insolation', np.sum( np.mean( insol, axis=1 ) * np.cos( np.deg2rad(lats[0,:]) ) ) / np.sum( np.cos( np.deg2rad( lats[0,:]))))
    #print('insol shape:', insol.shape)
    #print('insol [day1]:', insol[:,175])
    #--------- millenial scale with orbital params ---------
    #'''
    #   insol       = daily_insolation( lats[0,:], days, orb) # for having millenial scale insolation
    #   Qann    = np.mean(insol, axis=1)  # time average over the year
    #   #print('Qann.shape', Qann.shape)
    #   #print('Lat:', lats[0,:])
    #   Qglobal = np.empty_like( kyears )
    #   for n in range( kyears.size ):   # global area-weighted average
    #      Qglobal[n] = np.sum( Qann[:,n] * np.cos( np.deg2rad(lats[0,:]) ) ) / np.sum( np.cos( np.deg2rad(lats[0,:])))
    #      #print(Qglobal.shape)
    #   #print('rows and columns of lats',len(lats), len(lats[0]))
    #   #print('rows and columns of insol',len(insol), len(insol[0]))
    #   print('size of insol',insol.size)
    #   # To set up the model with different orbital parameters: orb = {'ecc':0., 'obliquity':90., 'long_peri':0.}
    #   orb_0  = OrbitalTable.interp(kyear=0)     # present-day orbital parameters
    #   orb_10 = OrbitalTable.interp(kyear=-10)   # orbital parameters for 10 kyrs before present
    #   orb_23 = OrbitalTable.interp(kyear=-23)   # 23 kyrs before present
    #   Q_0    = daily_insolation( lats[0,:], days, orb_0 )
    #   Q_10   = daily_insolation( lats[0,:], days, orb_10 )   # insolation arrays for each of the three sets of orbital parameters
    #   Q_23   = daily_insolation( lats[0,:], days, orb_23 )
    #   Qdiff  = Q_10 - Q_23
    #   print('The area-weighted global average of the difference:',np.average(np.mean(Qdiff,axis=1), weights=np.cos(np.deg2rad(lats[0,:]))))
    #'''
    #------------- *** End of insolation *** --------

    # Tensors which go into the state_vector
    global u1,h1,u2,h2,q1,b1,b2,q2,w1,w2
    u1   = sph.TensorField(1,S,domain) # wind speed lower layer [0]=meridional [1]=zonal
    h1   = sph.TensorField(0,S,domain)
    u2   = sph.TensorField(1,S,domain) # wind speed upper layer [0]=meridional [1]=zonal
    h2   = sph.TensorField(0,S,domain)
    q1   = sph.TensorField(0,S,domain)
    b1   = sph.TensorField(0,S,domain)
    b2   = sph.TensorField(0,S,domain)
    q2   = sph.TensorField(0,S,domain)
    w1   = sph.TensorField(0,S,domain)
    w2   = sph.TensorField(0,S,domain)

    global u10,u20,b10,b20,h10
    u10  = sph.TensorField(1,S,domain)
    u20  = sph.TensorField(1,S,domain)
    b10  = sph.TensorField(0,S,domain)
    b20  = sph.TensorField(0,S,domain)
    h10  = sph.TensorField(0,S,domain)

    global Du1,Du2,Db1,Db2,Dh1,Dh10,Dh2,uh1,uh2,uq1,uq2,uhf1,uqf1,uqf2,ubf1,uw1,uw2,uwf1,uwf2
    Du1  = sph.TensorField(2,S,domain)
    Du2  = sph.TensorField(2,S,domain)
    Db1  = sph.TensorField(1,S,domain)
    Db2  = sph.TensorField(1,S,domain)
    Dh1  = sph.TensorField(1,S,domain)
    Dh10 = sph.TensorField(1,S,domain)
    Dh2  = sph.TensorField(1,S,domain)
    uh1  = sph.TensorField(1,S,domain)
    uh2  = sph.TensorField(1,S,domain)
    uq1  = sph.TensorField(1,S,domain)
    uq2  = sph.TensorField(1,S,domain)
    #ub1  = sph.TensorField(1,S,domain)
    uhf1 = sph.TensorField(1,S,domain)
    uqf1 = sph.TensorField(1,S,domain)
    uqf2 = sph.TensorField(1,S,domain)
    ubf1 = sph.TensorField(1,S,domain)
    uw1  = sph.TensorField(1,S,domain)
    uw2  = sph.TensorField(1,S,domain)
    uwf1 = sph.TensorField(1,S,domain)
    uwf2 = sph.TensorField(1,S,domain)

    global divuh1,divuh2,divuq1,divuq2,divuw1,divuw2
    #global divub1,divuhf1,divuqf1,divuwf2,divubf1
    divuh1  = sph.TensorField(0,S,domain)
    divuh2  = sph.TensorField(0,S,domain)
    divuq1  = sph.TensorField(0,S,domain)
    divuq2  = sph.TensorField(0,S,domain)
    divuw1  = sph.TensorField(0,S,domain)
    divuw2  = sph.TensorField(0,S,domain)
    #divub1  = sph.TensorField(0,S,domain)
    #divuhf1 = sph.TensorField(0,S,domain)
    #divuqf1 = sph.TensorField(0,S,domain)
    #divuwf2 = sph.TensorField(0,S,domain)
    #divubf1 = sph.TensorField(0,S,domain)

    global u1_rhs,h1_rhs,u2_rhs,h2_rhs,q1_rhs,b1_rhs,b2_rhs,q2_rhs,w1_rhs,w2_rhs
    u1_rhs = sph.TensorField(1,S,domain)
    h1_rhs = sph.TensorField(0,S,domain)
    u2_rhs = sph.TensorField(1,S,domain)
    h2_rhs = sph.TensorField(0,S,domain)
    q1_rhs = sph.TensorField(0,S,domain)
    b1_rhs = sph.TensorField(0,S,domain)
    b2_rhs = sph.TensorField(0,S,domain)
    q2_rhs = sph.TensorField(0,S,domain)
    w1_rhs = sph.TensorField(0,S,domain)
    w2_rhs = sph.TensorField(0,S,domain)

    global Ev,Oc,CC1,DD1,DDr,VV2,VV1,Prec1,speed,R2,Rad1,Rad2,RT1,RT2,q_o_h1,q_o_h2
    Ev     = sph.TensorField(0,S,domain) # evaporation
    Oc     = sph.TensorField(0,S,domain) # Ocean
    CC1    = sph.TensorField(0,S,domain) # Condensation/Latent heat release forcing/ep
    DD1    = sph.TensorField(0,S,domain) # Downdraft+bias correction due to condensation in the lower layer, from upper troposphere to the lower troposphere
    DDr    = sph.TensorField(0,S,domain) # Downdraft without bias correction due to condensation in the lower layer, from upper troposphere to the lower troposphere
    VV2    = sph.TensorField(0,S,domain) # Vaporization in the upper layer
    VV1    = sph.TensorField(0,S,domain) # Vaporization in the lower layer
    Prec1  = sph.TensorField(0,S,domain) # Precipitation in the lower layer
    speed  = sph.TensorField(0,S,domain) # velocity magnitude
    R2     = sph.TensorField(1,S,domain)
    Rad1   = sph.TensorField(0,S,domain) # Newtonian cooling for the lower layer
    Rad2   = sph.TensorField(0,S,domain) # Newtonian cooling for the upper layer
    RT1    = sph.TensorField(0,S,domain) # RRTM_LW/SW, lower layer
    RT2    = sph.TensorField(0,S,domain) # RRTM_LW/SW, upper layer
    q_o_h1 = sph.TensorField(0,S,domain) # q1/h1
    q_o_h2 = sph.TensorField(0,S,domain)

    global Press1, Press2
    Press1 = np.empty((ni,nj),dtype=np.float64)
    Press2 = np.empty((ni,nj),dtype=np.float64)

    global state_vector,RHS,timestepper
    state_vector = eq.StateVector(u1,h1,u2,h2,q1,b1,b2,q2,w1,w2)
    RHS          = eq.StateVector(u1,h1,u2,h2,q1,b1,b2,q2,w1,w2)
    timestepper  = timesteppers.SBDF2(eq.StateVector, u1,h1,u2,h2,q1,b1,b2,q2,w1,w2)

    global t
    t = 0

    #
    # now read initial conditions
    #
    # topography height and ocean mask were already passed as parameters,
    # matching the selected grid resolution/size
    h10['g'][0] = (delta1)*itopo*HORO/(gamma_topo*g_real*H0_dim)
    Oc['g'][0]  = 1-land_frac

    #
    # if there is a restart file, read it.
    #
    restart_file = None
    if os.path.exists("INPUT/aeolus2.res.nc"):
        print("reading restart file INPUT/aeolus2.res.nc")
        restart_file = Dataset("INPUT/aeolus2.res.nc")
    elif os.path.exists("INPUT/aeolus2_res.npz"):
        print("reading restart file INPUT/aeolus2_res.npz")
        restart_file = np.load("INPUT/aeolus2_res.npz")

    # Reading variables and cutting out the local domain should work
    # for both file formats in the same way.
    if (restart_file):
        init_warmstart(restart_file)
    else:
        init_coldstart()

    # build matrices
    global P,M,L,LU
    m_start = domain.distributor.coeff_layout.start(1)[0]
    m_len   = domain.distributor.coeff_layout.local_shape(1)[0]
    m_end   = m_start + m_len - 1
    for m in range(m_start,m_end+1):
        Mm,Lm = eq.shallow_water(S,m,[g0,B3,B2,B1,H1,H2,Om,a,nu,kappa,Q01,Q02])
        M.append(Mm.astype(np.complex128))
        L.append(Lm.astype(np.complex128))
        P.append(0.*Mm.astype(np.complex128))
        #LU.append([None])

    #
    # prepare diagnostic output
    #
    # TODO convert outpout variables to real-world units
    # Output variables are nondimensional, the corresponding dimensional scales are mentioned below. Dimensional scales may vary with respect to nondimensionalization method.
    # dimension of velocity on the equatorial beta plane is beta*L_d^2, where L_d is equatorial Rossby deformation radius
    # layer 1: lower layer, layer 2: upper layer
    # q presets bulk of humidity and its scale is [L_v*g*H1/(C_p*theta_s*45)], however it can be rescaled to vertically averaged specific humidity with [(1/45)Kg/Kg] scale.
    # Time scale is T=1/(beta*L_d^2) if x,y are nondimentionalized by equatorial Rossby deformation radius.
    # RT: Thermal forcing due to Radiative Transfer
    # define and write out topography as static 2-D variable
    ncout.AddVariable('h10', 'effective topography height',                  '[gamma_topo*g_real*H0_dim]',               miss_val=None, axes=axes_yx)
    ncout.AddVariable('ocean', 'land-sea mask',                              '1',                                        miss_val=None, axes=axes_yx)
    ncout.AddVariable('u1th', 'Meridional velocity layer 1',                 '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('u1ph', 'Zonal (azimuthal) velocity layer 1',          '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('h1',   'Pseudo-height layer 1',                       'H',                                        miss_val=None, axes=axes_tyx)
    ncout.AddVariable('q1',   'Bulk of Specific humidity at layer 1',        '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
    ncout.AddVariable('q2',   'Bulk of Specific humidity at layer 2',        '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
    ncout.AddVariable('w1',   'Bulk of Precipitable Water at layer 1',       '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
    ncout.AddVariable('w2',   'Bulk of Precipitable Water at layer 2',       '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
    ncout.AddVariable('Prec1','Bulk of Precipitaion at layer 1',             '[(L_v.g/(C_p.theta_s))Kg/Kg]',             miss_val=None, axes=axes_tyx)
    ncout.AddVariable('b1',   'Buoyancy layer 1',                            '[g*theta/theta_s]',                        miss_val=None, axes=axes_tyx)
    ncout.AddVariable('b2',   'Buoyancy layer 2',                            '[g*theta/theta_s]',                        miss_val=None, axes=axes_tyx)
    ncout.AddVariable('CC1',  'CLWC: Cloud Liquid Water Content at layer 1', '[Q/T]',                                    miss_val=None, axes=axes_tyx)
    ncout.AddVariable('DD1',  'Downdraft to layer 1 (balanced)',             '[Q/T]',                                    miss_val=None, axes=axes_tyx)
    ncout.AddVariable('DDr',  'Downdraft to layer 1 (unbalanced)',           '[Q/T]',                                    miss_val=None, axes=axes_tyx)
    ncout.AddVariable('Ev',   'Sea surface evaoporation',                    '[Q/T]',                                    miss_val=None, axes=axes_tyx)
    ncout.AddVariable('u2th', 'Meridional velocity layer 2',                 '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('u2ph', 'Zonal (azimuthal) velocity layer 2',          '[(g*H)^0.5], [beta*L_d^2] at the equator', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('h2',   'Pseudo-height layer 2',                       'H',                                        miss_val=None, axes=axes_tyx)
    ncout.AddVariable('RT1',   'Radiative transfer flux, lower layer',       '[HB]',                                     miss_val=None, axes=axes_tyx)
    ncout.AddVariable('RT2',   'Radiative transfer flux, upper layer',       '[HB]',                                     miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('R2th', 'long name of R2th', 'unknown units', miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('R2ph', 'long name of R2ph', 'unknown units', miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('divuhf1', 'long name of divuhf1', 'unknown units', miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('divuqf1', 'long name of divuqf1', 'unknown units', miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('divuwf1', 'long name of divuwf1', 'unknown units', miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('divuwf2', 'long name of divuwf2', 'unknown units', miss_val=None, axes=axes_tyx)
    #ncout.AddVariable('divubf1', 'long name of divubf1', 'unknown units', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('menthalpy', 'menthalpy = h1+H1-ep1*q1-ep1*Q01', 'unknown units', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('insolation', 'daily mean insolation', 'W/m2?', miss_val=None, axes=axes_ty)
    ncout.AddVariable('Press1', 'Pressure at lower layer (P_level1)', '[Pa]', miss_val=None, axes=axes_tyx)
    ncout.AddVariable('Press2', 'Pressure at upper layer (P_level2)', '[Pa]', miss_val=None, axes=axes_tyx)


    # topography and land-sea mask were already written out as static variables above.
    #ncout.WriteVar('h10',  h10['g'][0],    t*T_scale)
    #ncout.WriteVar('ocean', Oc['g'][0],    t*T_scale)

    # enable strict floating point exceptions to catch errors/instabilities as early as possible
    np.seterr(all='raise')
    print("aeolus2_init() finished", flush=True)

# end of aeolus2_init()


def init_coldstart():
    """
    TODO: Initialise with synthetic/idealized conditions.
    And/or read idealized initial conditions from climatology files.
    """
    global u1,h1,u2,h2,q1,b1,b2,q2,w1,w2
    #global u10,u20,b10,b20,h10

    # TODO: put preprocessing code into here.
    # The parameters used in preprocessing must match the parameters
    # used at model run time. If preprocessing is done in a separate script,
    # that opens up an entire universe of possibilities for inconsistencies,
    # which lead to strange model behaviour without any hint the reasons.
    #
    # Input: a (NetCDF) file with 3-D values for u-wind, v-wind and temperature,
    # on pressure levels,
    # aggregated to one time step (daily mean, monthly mean, whatever),
    # and interpolated to the model grid.
    # cdo remapbil\,POEM/exp/MOM_LAD_AEOLUS2/griddes-Aeolus2-192x96.txt \
    #     -selday\,1 \
    #     -ydaymean \
    #     -cat \
    #     -unpack \
    #     -copy '198*/uvt_EI_daily_mean_198*_01.nc' \
    #    eraint-uvt-1-1000-ydaymean-01-jan-1980-1989.nc
    #
    # cdo remapbil\,POEM/exp/MOM_LAD_AEOLUS2/griddes-Aeolus2-192x96.txt \
    #     -monmean \
    #     -unpack \
    #     -sellevel\,1000\,900\,800\,700\,500\,400\,300 \
    #      '1980/uvt_EI_daily_mean_1980_01.nc' \
    #    eraint-uvt-7levels-1000-monmean-jan-1980.nc
    #
    #
    # Parameter delta1 specifies the aspect ratio, i.e. split of
    # pressure levels onto the 2 model layers.
    # Troposphere ranges from ca. 1000 mbar to ca. 100 mbar.
    # Thus: count number of levels in input file which are below or equal to (i.e. >= ) 100 mbar.
    # Split that number according to delta1.
    # Compute mean of lower and upper layer pressure values as representing numbers for the model layers,
    # p_level1 and p_level2.

    # Look into TWO_L_ERIM_to_Aeolus2.m for the processing steps.

    ## 0th step: read from the input file the data slices for the
    # local data domain [ics:ice,jcs:jce]

    ## 0.0 open the input file
    # TODO: make the name of the input file a runtime parameter
    #infilename = 'INPUT/eraint-uvt-1-1000-ydaymean-01-jan-1980-1989.nc'
    infilename = 'INPUT/eraint-uvt-7levels-monmean-jan-1980.nc'
    infile = Dataset(infilename, 'r', parallel=False)
    print("init_coldstart() Reading ",infilename," for initial conditions", flush=True)

    ## 0.1 identify lon,lat,level,time dimensions
    for v in infile.variables:
        if v in ('longitude','lon','X','x'):
            lonaxisname = v
            lonaxisvar = infile.variables[v]
            print("found longitude axis: ", v)
        if v in ('latitude','lat','Y','y'):
            lataxisname = v
            lataxisvar = infile.variables[v]
            print("found latitude axis: ", v)
        if v in ('z','Z','level','Level','LEVEL'):
            levaxisname = v
            levaxisvar = infile.variables[v]
            print("found level axis: ", v)
        if v in ('time','Time','TIME'):
            #timaxisname = v
            timaxisvar = infile.variables[v]
            print("found time axis: ", v)
        if v in ('t','T','temp','T3D'):
            #tvarname = v
            tvar = infile.variables[v]
            print("found temperature: ", v)
        if v in ('u','U','uwind','Uwind'):
            #uvarname = v
            uvar = infile.variables[v]
            print("found u-wind: ", v)
        if v in ('v','V','vwind','Vwind'):
            #vvarname = v
            vvar = infile.variables[v]
            print("found v-wind: ", v)

    ## 0.2 check that size and axis orientation of the horizontal grid of the input
    #     file matches the model grid
    if Nxx != lonaxisvar.size:
        print("length of longitude axis of initial conditions files does not match model grid size: ",
              Nxx, " != ", lonaxisvar.size)
        raise Exception
    if not lonaxisvar.units in ("degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE"):
        print("unknown unit ",lonaxisvar.units," for longitude axis variable",lonaxisname)
        raise Exception
    if not np.any(double_eq(lonaxisvar[ics:ice],lons_local[:])):
        print("Error: longitude axis of initial conditions file does not match model grid")
        print("from initial conditions:",lonaxisvar[ics:ice])
        print("from model grid:",lons_local[:])
        #print("match:",lonaxisvar[:] != lons_local[ics:ice])
        #print("double_eq:",double_eq(lonaxisvar[ics:ice],lons_local[:]))
        print("diff:",lonaxisvar[ics:ice]-lons_local[:])
        raise Exception
    if Nyy != lataxisvar.size:
        print("length of latitude axis of initial conditions files does not match model grid size: ",
              Nyy, " != ", lataxisvar.size)
        raise Exception
    if not lataxisvar.units in ("degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN"):
        print("unknown unit ",lataxisvar.units," for latitude axis variable",lataxisname)
        raise Exception
    if not np.any(double_eq(lataxisvar[jcs:jce],lats_local[:],1e-5)):
        # Wow. really large rounding errors might occur here. But what can we do...
        print("Error: latitude axis of initial conditions file does not match model grid")
        print("from initial conditions:",lataxisvar[jcs:jce])
        print("from model grid:",lats_local[:])
        print("diff:",lataxisvar[jcs:jce]-lats_local[:])
        raise Exception

    if timaxisvar is None:
        print("No time axis found")
    else:
        if timaxisvar.size > 1:
            print("Warning: input file has more than one timestep.")
            print("         will use only the first one for initial conditions.")

    ## 0.3 check direction and units of vertical axis
    #     [mbar], [millibars], [hPa], [Pa]; upward or downward
    #     Convert values to [Pa]
    if not levaxisvar.units in ("mbar","millibars","hPa","Pa"):
        print("unknown unit ",levaxisvar.units," for level axis variable",levaxisname)
        raise Exception
    # SIGH. range() & Co produce indices one beyond the actual array
    # size.  Which might be convenient for people who do not
    # understand for-loops in C, But is a pain in the ass for
    # reversing the iteration direction in a generic way.  (well, or
    # maybe I am just too stupid to appreciate this feature)
    if levaxisvar[0] < levaxisvar[-1]:
        print("level axis downwards from ", levaxisvar[0], " to ", levaxisvar[-1])
        indbot = levaxisvar.size-1; indtop = -1; levstep = -1
    else:
        print("level axis upwards from ", levaxisvar[0], " to ", levaxisvar[-1])
        indbot = 0; indtop = levaxisvar.size; levstep = 1

    ## 0.4 check units of input variables
    #     [m/s] for u,v wind, [K] or [C] for t
    if not tvar.units in ("K", "degK", "C", "degC", "deg C"):
        print("unknown unit ",tvar.units," for temperature variable",tvar,
              ". Need Kelvin [K] or Celsius [C]")
        raise Exception
    if not uvar.units in ("m/s", "m s**-1"):
        print("unknown unit ",uvar.units," for u-wind variable",uvar)
        raise Exception
    if not vvar.units in ("m/s", "m s**-1"):
        print("unknown unit ",vvar.units," for v-wind variable",vvar)
        raise Exception


    ## 1st step: vertical integration (average) from pressure levels onto model layers
    # We must copy the level axis, so that we have write permissions
    # in case we need to convert units.

    ## 1.0 find limits and border indices of input levels for the model layers.
    levunits = levaxisvar.units
    levvar = levaxisvar[:] # Pa
    if levunits in ("mbar","millibars","hPa"):
        levvar[:] = levaxisvar[:] * 100 # convert to Pa
        levunits = "Pa"
    print("delta1 ", delta1, " split between lower and upper model layer") # from config
    # Sigh. delta1 = 0.6 means 60% of troposphere should go to lower model layer.
    # But since pressure levels have lower numbers for higher levels, we need to revert.
    # Also be careful with rounding errors here.
    borderlevel = 100000.0*(1.0-delta1) # [Pa]
    print("-> lower model layer between 100000 Pa and ", borderlevel, " Pa")
    ntropolevels = 0
    nlowerlevels = 0
    for i in range(indbot, indtop, levstep):
        print(i, levaxisvar[i], levvar[i], borderlevel)
        if levvar[i] >= 10000: # Pa
            ntropolevels += 1
            indtopupper = i
        if levvar[i] > borderlevel:
            nlowerlevels += 1
            indtoplower = i
    print("found ", ntropolevels, " levels for the troposphere")
    print("found ", nlowerlevels, " levels for the lower model layer")
    print("lower layer from ", levvar[indbot], " to ", levvar[indtoplower], " Pa")
    nupperlevels = ntropolevels - nlowerlevels
    indbotupper = indtoplower+levstep
    print("found ", nupperlevels, " levels for the upper model layer")
    print("upper layer from ", levvar[indbotupper], " to ", levvar[indtopupper], " Pa")

    # SIGH again.
    # More printout-loops to make sure that we got all the bottom and top indices right.
    #print("input levels for the lower model layer:")
    #for i in range(indbot,indtoplower+levstep,levstep):
    #    print(i, levaxisvar[i])
    #print("input levels for the upper model layer:")
    #for i in range(indbotupper,indtopupper+levstep,levstep):
    #    print(i, levaxisvar[i])

    ## 1.1 the actual vertical integration of wind speeds
    #print("shape uvar", uvar.shape)
    # this should yield 3 dimensions (lev,lat,lon) of apropriate sizes
    # TODO: test that, and raise exeception in case of mismatch
    if timaxisvar is None:
        print("shape of uvar for lower layer:",
              uvar[indbot:indtoplower:levstep,jcs:jce,ics:ice].shape)
        print("shape of uvar for upper layer:",
              uvar[indbotupper:indtopupper:levstep,jcs:jce,ics:ice].shape)
        ulow = np.mean(uvar[indbot:indtoplower:levstep,jcs:jce,ics:ice], axis=0)
        uupp = np.mean(uvar[indbotupper:indtopupper:levstep,jcs:jce,ics:ice], axis=0)
        vlow = np.mean(vvar[indbot:indtoplower:levstep,jcs:jce,ics:ice], axis=0)
        vupp = np.mean(vvar[indbotupper:indtopupper:levstep,jcs:jce,ics:ice], axis=0)
    else:
        print("shape of uvar for lower layer:",
              uvar[0,indbot:indtoplower:levstep,jcs:jce,ics:ice].shape)
        print("shape of uvar for upper layer:",
              uvar[0,indbotupper:indtopupper:levstep,jcs:jce,ics:ice].shape)
        ulow = np.mean(uvar[0,indbot:indtoplower:levstep,jcs:jce,ics:ice], axis=0)
        uupp = np.mean(uvar[0,indbotupper:indtopupper:levstep,jcs:jce,ics:ice], axis=0)
        vlow = np.mean(vvar[0,indbot:indtoplower:levstep,jcs:jce,ics:ice], axis=0)
        vupp = np.mean(vvar[0,indbotupper:indtopupper:levstep,jcs:jce,ics:ice], axis=0)

    ## 1.2 convert temperature to potential temperature

    # In TWO_L_ERIM_to_Aeolus2.m we find:
    #% R/Cp for conversion of t to potential temperature
    #R_o_Cp          = 0.285;
    #% convert temperature to potential temperature
    #T_2D_300_i      = squeeze(mean(T_4D(:,:,1,ti:tf),4))*((1000/300) .^R_o_Cp);
    # [where mean(X(...,ti:tf),4) is mean over time]
    # (Note the scaling of pressure heights!)

    # From
    # https://en.wikipedia.org/wiki/Potential_temperature
    # The potential temperature of a parcel of fluid at pressure $P$
    # is the temperature that the parcel would attain if adiabatically
    # brought to a standard reference pressure $P_0$, usually 1,000
    # hPa (1,000 mb). The potential temperature is denoted
    # $\theta=T\left({\frac {P_{0}}{P}}\right)^{R/c_{p}}$ , where $T$$
    # is the current absolute temperature (in K) of the parcel, $R$ is
    # the gas constant of air, and $c_p$ is the specific heat capacity
    # at a constant pressure. $R/c_p = 0.286$ for air (meteorology).
    #
    # Sigh. The name theta is in Dedalus already used for (co-)latitude
    # coordinates.

    # MOM5 / FMS constants.F90 says:
    #! <DATA NAME="RDGAS" UNITS="J/kg/deg" TYPE="real" DEFAULT="287.04">
    #!   gas constant for dry air
    #! </DATA>
    #! <DATA NAME="KAPPA" TYPE="real" DEFAULT="2./7.">
    #!   RDGAS / CP_AIR
    #! </DATA>
    #! <DATA NAME="CP_AIR" UNITS="J/kg/deg" TYPE="real" DEFAULT="RDGAS/KAPPA">
    #!   specific heat capacity of dry air at constant pressure
    #! </DATA>
    # Thus:
    # RDGAS/CP_AIR = RDGAS/(RDGAS/KAPPA) = KAPPA = 2/7 = 0.285714

    kappa = 2.0/7.0
    if tvar.units in ("C", "degC", "deg C"):
        tcorr = 273.15
    else:
        tcorr = 0.0

    tpot = np.empty((levaxisvar.size, nj, ni), dtype=np.float64)
    for i in range(indbot, indtop, levstep):
        if timaxisvar is None:
            tpot[i,:,:] = np.power((tvar[i,jcs:jce,ics:ice]+tcorr)*(100000.0/levvar[i]), kappa)
        else:
            tpot[i,:,:] = np.power((tvar[0,i,jcs:jce,ics:ice]+tcorr)*(100000.0/levvar[i]), kappa)
    print("shape of tpot ", tpot.shape)

    ## 1.3 the actual vertical integration of potential temperature
    tpotlow = np.mean(tpot[indbot:indtoplower:levstep,:,:], axis=0)
    tpotupp = np.mean(tpot[indbotupper:indtopupper:levstep,:,:], axis=0)
    print("shape of tpotlow, tpotupp ", tpotlow.shape, tpotupp.shape)


    ## 2nd step: convert u,v,tpot on model layers to non-dimensional quantities

    ## 2.1 scaling wind speeds
    #g                 = 9.8;
    #H0                = 10000;%17500;%20000;% https://www.e3s-conferences.org/articles/e3sconf/pdf/2019/02/e3sconf_icst2018_04002.pdf
    #%H0              = 28000;% if we choose Hoskin pseudo-height
    #% Velocity scale for converting wind speeds
    #U_bt              = sqrt(g*H0);                % 313ms^-1, also U scale
    #U_scale         = U_bt;
    print("U_scale ", U_scale, " for wind speeds scaling") # from config file

    #% non-dimensional velocities for lower and upper model layer (on 2-D maps)
    #
    #U_2D_low  = (U_2D_500_i+U_2D_700_i+U_2D_800_i+U_2D_900_i+U_2D_1000_i)/5/U_bt;
    #U_2D_upp  = (U_2D_300_i+U_2D_400_i)/2/U_bt;
    #V_2D_low  = (V_2D_500_i+V_2D_700_i+V_2D_800_i+V_2D_900_i+V_2D_1000_i)/5/U_bt;
    #V_2D_upp  = (V_2D_300_i+V_2D_400_i)/2/U_bt;
    u1['g'][1] = np.transpose(ulow) / U_scale
    u2['g'][1] = np.transpose(uupp) / U_scale
    u1['g'][0] = np.transpose(vlow) / U_scale
    u2['g'][0] = np.transpose(vupp) / U_scale

    ## 2.2 scaling potential temperature
    # T_s0 is derived as min(min(T_2D_low_D)).
    # *** to compare comparables B2, TS, should be the same for all months **
    # 01-Jan-1980 : 250.0277, monthly mean Jan 1980 : 262.1009 ERA-Interim
    # T_s0          = 250.0277;
    print("T_s0 ", T_s0, " for potential temperature scaling") # from config file
    # deltab      = (mean(T_2D_upp_D(:))- mean(T_2D_low_D(:)))/T_s0;
    #deltab = (np.mean(tpotupp) - np.mean(tpotlow)) / T_s0
    B1 = np.mean(tpotlow) / T_s0
    B2 = np.mean(tpotupp) / T_s0
    print("potential temperature mean for upper and lower layer ",
          np.mean(tpotupp), " ", np.mean(tpotlow))
    print("scaled with T_s0: B1 ", B1, " B2 ", B2)
    #print("deltab ", deltab, " for potential temperature scaling")
    #% convert potential temperature to non-dimensional temperature anomalies
    #T_2D_low    = (T_2D_low_D-TS)/TS;
    #T_2D_upp    = ((T_2D_upp_D-TS)/TS)- deltab;
    #T_1D_low    = mean(T_2D_low);
    #T_1D_upp    = mean(T_2D_upp);
    #B2          = 1.1514;% mean(T_2D_upp_D(:))/mean(T_2D_low_D(:))+deltab;
    #T_2D_low    = T_2D_low_D/TS - B1;
    #T_2D_upp    = T_2D_upp_D/TS - B2;%- deltab;
    b1['g'][0][:,:] = np.transpose(tpotlow)/T_s0 - B1
    b2['g'][0][:,:] = np.transpose(tpotupp)/T_s0 - B2 # - deltab

    ## 2.3 write the preprocessed data into a diagnostic NetCDF file
    outfile = NetCDFOutput(comm,
                           coords = {
                               "nlons_global" : Nxx,
                               "nlats_global" : Nyy,
                               "lons"         : lons_local[:],
                               #"lonb"         : lonb[:],
                               "lats"         : lats_local[:],
                               #"latb"         : latb[:],
                               "ics"          : ics,
                               "ice"          : ice, # according to python conventions, one after end
                               "jcs"          : jcs,
                               "jce"          : jce, # according to python conventions, one after end
                           },
                           f = os.path.join(output_folder, "a2coldstartdata.nc"),
                           comment = "Results of preprocessing coldstart data derived from "
                           + infilename)
    axes_yx  = (ncout.lataxisname,ncout.lonaxisname,)  # static maps
    outfile.AddVariable('U_Scale', 'for wind speed scaling', '1', miss_val=None, axes=())
    outfile.AddVariable('T_s0', 'for potential temperature scaling', '1', miss_val=None, axes=())
    #outfile.AddVariable('deltab', 'for potential temperature scaling', '1', miss_val=None, axes=())
    # We better avoid variables which are distinguished only by upper vs lower case names,
    # to prevent confusion in tools like ferret.
    # Thus append suffix _scale to these names.
    outfile.AddVariable('B1_scale', 'mean potential temperature on lower layer scaled with T_s0', '1', miss_val=None, axes=())
    outfile.AddVariable('B2_scale', 'mean potential temperature on upper layer scaled with T_s0', '1', miss_val=None, axes=())
    outfile.AddVariable('U1', 'zonal wind at lower model layer scaled nondimensionally', '1', miss_val=None, axes=axes_yx)
    outfile.AddVariable('U2', 'zonal wind at upper model layer scaled nondimensionally', '1', miss_val=None, axes=axes_yx)
    outfile.AddVariable('V1', 'meridional wind at lower model layer scaled nondimensionally', '1', miss_val=None, axes=axes_yx)
    outfile.AddVariable('V2', 'meridional wind at upper model layer scaled nondimensionally', '1', miss_val=None, axes=axes_yx)
    outfile.AddVariable('b1', 'potential temperature at lower model layer scaled nondimensionally', '1', miss_val=None, axes=axes_yx)
    outfile.AddVariable('b2', 'potential temperature at upper model layer scaled nondimensionally', '1', miss_val=None, axes=axes_yx)

    outfile.WriteVar('U_Scale', U_scale, timestamp=None)
    outfile.WriteVar('T_s0', T_s0, timestamp=None)
    #outfile.WriteVar('deltab', deltab, timestamp=None)
    outfile.WriteVar('B1_scale', B1, timestamp=None)
    outfile.WriteVar('B2_scale', B2, timestamp=None)
    outfile.WriteVar('U1', u1['g'][1], timestamp=None)
    outfile.WriteVar('U2', u2['g'][1], timestamp=None)
    outfile.WriteVar('V1', u1['g'][0], timestamp=None)
    outfile.WriteVar('V2', u2['g'][0], timestamp=None)
    outfile.WriteVar('b1', b1['g'][0], timestamp=None)
    outfile.WriteVar('b2', b2['g'][0], timestamp=None)

    outfile.close()

    ## 42nd step init the remaining variables

    h1['g'][0][:,:]     = 0.0
    h2['g'][0][:,:]     = 0.0
    #u1['g'][0][:,:]     = 0.0 # lower layer v-wind
    #u1['g'][1][:,:]     = 0.0 # lower layer u-wind
    #u2['g'][0][:,:]     = 0.0 # upper layer u-wind
    #u2['g'][1][:,:]     = 0.0 # upper layer v-wind
    q1['g'][0][:,:]     = 0.0
    #b1['g'][0][:,:]     = 0.0
    #b2['g'][0][:,:]     = 0.0
    q2['g'][0][:,:]     = 0.0
    b10['g'][0][:,:]     = 0.0
    b20['g'][0][:,:]     = 0.0
    #w1['g'][0][:,:]     = 0.0
    #w2['g'][0][:,:]     = 0.0

# end of init_coldstart()


def init_warmstart(restart_file):
    """
    TODO plausibility checks: grid size etc
    Not todo here: scaling, application of forcing, etc.
    A run with save/restart should yield exactly
    the same results as a run without intermediate restart.
    """
    global u1,h1,u2,h2,q1,b1,b2,u10,u20,b10,b20,h10,q2,w1,w2,t

    #print(restart_file['h1'])
    print('h1 shape ', restart_file['h1'].shape)
    h1['g'][0]          = restart_file['h1'][ics:ice,jcs:jce]
    h2['g'][0]          = restart_file['h2'][ics:ice,jcs:jce]
    u1['g'][0]          = restart_file['u1th'][ics:ice,jcs:jce]
    u1['g'][1]          = restart_file['u1ph'][ics:ice,jcs:jce]
    u2['g'][0]          = restart_file['u2th'][ics:ice,jcs:jce]
    u2['g'][1]          = restart_file['u2ph'][ics:ice,jcs:jce]
    q1['g'][0]          = restart_file['q1'][ics:ice,jcs:jce]
    b1['g'][0]          = restart_file['b1'][ics:ice,jcs:jce]
    b2['g'][0]          = restart_file['b2'][ics:ice,jcs:jce]
    q2['g'][0]          = restart_file['q2'][ics:ice,jcs:jce]
    b10['g'][0]         = restart_file['b1'][ics:ice,jcs:jce]
    b20['g'][0]         = restart_file['b2'][ics:ice,jcs:jce]
    #w1['g'][0]          = 0.*restart_file['q1'][ics:ice,jcs:jce]#restart_file['w1'][ics:ice,jcs:jce]
    #w2['g'][0]          = 0.*restart_file['q1'][ics:ice,jcs:jce]#restart_file['w2'][ics:ice,jcs:jce]
    t                   = restart_file['t'][0]
# end of init_warmstart()


def aeolus2_finish(Time_seconds, Time_days):
    """
    write out restart files and clean up
    """
    aeolus2_restart(Time_seconds, Time_days)
    # TODO: implement aeolus2_finish()
# end of aeolus2_finish()


def aeolus2_restart(Time_seconds, Time_days):
    """
    write out (intermediate) restart files
    """
    # TODO: implement aeolus2_restart()
    # Write out a restart file that can be read back
    # in function init_warmstart()
# end of aeolus2_restart()


ISTOCK_WATER=1
ISTOCK_HEAT=2
ISTOCK_SALT=3
def aeolus2_get_stock_pe(idx):
    """
    Compute current stock of water/heat/salt present in the local data
    domain of this processing element.
    If checking of mass/energy conservation is switched on in the model
    configuration, this is called once after initialisation, and
    then after each time step.
    idx: index specifying which stock to compute
    return value: scalar output parameter
    """
    # TODO: Masoud+Sullyandro: implement the stock taking computations
    if idx == ISTOCK_WATER:
        # TODO: implement local water stock taking [kg]
        print("TODO: implement local water stock taking")
        return 0.0
    elif idx == ISTOCK_HEAT:
        # TODO: implement local heat stock taking [J]
        print("TODO: implement local heat stock taking")
        return 0.0
    elif idx == ISTOCK_SALT:
        # no salt in the atmosphere
        return 0.0
    else:
        # unknown stock type
        return 0.0

# end of aeolus2_get_stock_pe()


def aeolus2_get_stock_atm_global(idx):
    """
    Compute current stock of water/heat/salt present globally.
    If checking of mass/energy conservation is switched on in the model.
    idx: index specifying which stock to compute
    return value: scalar output parameter
    """
    val = aeolus2_get_stock_pe(idx)
    return comm.allreduce(val, MPI.SUM)
# end of aeolus2_get_stock_atm_global()


def aeolus2_get_bottom_mass():
    """
    Return quantities at bottom of atmospere, as 2-D arrays.
    If CO2 is configured as tracer at the coupler interface,
    that value must also be returned.
    t_bot  near surface temperature [K]
    q_bot  near surface mixing ratio [kg/kg]
    p_bot  pressure at which atmos near surface values are assumed to be defined [pa]
    z_bot  height at which atmos near surface values are assumed to be defined [m]
           height above the surface for the lowest model level (m)
    p_surf surface pressure [pa]
    slp    sea level pressure [pa]
    co2_bot [kg/(kg wet air)]
    """
    t_bot   = np.empty((ni,nj),dtype=np.float64,order='F') # near surface temperature [K]
    q_bot   = np.empty((ni,nj),dtype=np.float64,order='F') # near surface mixing ratio [kg/kg]
    p_bot   = np.empty((ni,nj),dtype=np.float64,order='F') # pressure at which atmos near surface values are assumed to be defined [pa]
    # height above the surface for the lowest model level (m)
    z_bot   = np.empty((ni,nj),dtype=np.float64,order='F') # height at which atmos near surface values are assumed to be defined [m]
    p_surf  = np.empty((ni,nj),dtype=np.float64,order='F') # surface pressure [pa]
    slp     = np.empty((ni,nj),dtype=np.float64,order='F') # sea level pressure [pa]
    co2_bot = np.empty((ni,nj),dtype=np.float64,order='F') # [kg/(kg wet air)]
    # TODO: Masoud+Sullyandro: implement aeolus2_get_bottom_mass
    print("TODO: Masoud+Sullyandro: implement aeolus2_get_bottom_mass")

    # mass fraction of water vapour (relative humidity) at lower model layer
    f_mf = q1['g'][0]/(global_amax(comm, q1['g'][0])+0.0000001)
    # effective gas constant
    R_real = R_dry*(1.0 - f_mf) + R_moist*f_mf

    t_sl = T_s0*(b1['g'][0]+B1) - (g_real/R_real)*ma.log((1.0-delta1)/2.0) # Atmosphere temperature at sea level
    # adjust temperature at surface with respect to HORO
    gamma_trop = 0.00649 #  lapse rate of the troposphere [degree K/m]
    t_bot[:,:] = t_sl - gamma_trop * HORO

    asr_rh     = 1.0 # Aspect ratio for calibrating bulk relative humidity to its corresponding surface level.
    # TODO: adjust humidity with respect to HORO
    q_bot[:,:] = asr_rh*f_mf

    # pressure anomaly
    p_anom = g_real* ((h1['g'][0]+h2['g'][0]+h10['g'][0])*b1['g'][0])/a
    # pressure at lower model level
    # delta1 = 0.6 means lower layer from 1000 hPa to 400 hPa,
    # i.e. centered at / represented by 1000-(1000-400)/2
    Press1[:,:] = 1000*100 * (1+p_anom) # p_level1 == lowest pressure level 1000mb convert to Pa
    p_bot[:,:] = Press1

    #z_bot[:,:] = HORO # TODO: just for testing the communication
    z_bot[:,:] = 10.0 # copied from Aeolus 1.0

    p_surf[:,:] = Press1

    # TODO: adjust pressure to sealevel with respect to HORO
    slp[:,:] = 1000*100 * (1+p_anom) # P_s == lowest pressure level 1000mb convert to Pa

    if (co2_ind < 0):
        co2_bot[:,:] = 0.0
    # else:
    #   compute return co2 value for co2 cycle
    return t_bot, q_bot, p_bot, z_bot, p_surf, slp, co2_bot
# end of aeolus2_get_bottom_mass()


def aeolus2_get_bottom_wind():
    """
    returns u and v on the mass grid at the lowest model level
    u_bot [m/s]
    v_bot [m/s]
    """

    u_bot = np.empty((ni,nj),dtype=np.float64,order='F') # [m/s]
    v_bot = np.empty((ni,nj),dtype=np.float64,order='F') # [m/s]
    # TODO: Masoud+Sullyandro: implement aeolus2_get_bottom_wind
    print("TODO: Masoud+Sullyandro: implement aeolus2_get_bottom_wind", flush=True)
    # 15.2.24: Stefan just guessed a formula, just to have something to test the FMS interface.
    u_bot[:,:] = u1['g'][1]*U_scale
    v_bot[:,:] = u1['g'][0]*U_scale
    return u_bot, v_bot
# end of aeolus2_get_bottom_wind()


def aeolus2_update_down(x):
    """
    const int Time_seconds, Time_days, dayoftheyear

    Most of the following parameters are 2-D arrays,
    - except the tracer arrays denoted as 3D arrays -
    defined on the local data domain with index range (ics:ice,jcs:jce)

    Input parameters
    Quantities going from land+ice to atmos
    const double land_frac,       // fraction amount of land in a grid box [1]
    const double t_surf,          // surface temperature for radiation calculations [K]
    const double albedo,          // surface albedo for radiation calculations [1]
    const double albedo_vis_dir,
    const double albedo_nir_dir,
    const double albedo_vis_dif,
    const double albedo_nir_dif,
    const double rough_mom,       // surface roughness (used for momentum) [m]

    #DC: following block of input parameters are not required in Aeolus:
    #surface winds are directly calculated from roughness_length only.

    In the FMS coupler, subroutine sfc_boundary_layer(),
    u_star and b_star are defined so that u_star**2 = magnitude
    of surface stress divided by density of air at the surface,
    and u_star*b_star = buoyancy flux at the surface.

    const double u_star,          // friction velocity [m/s]
    const double b_star,          // bouyancy scale [m/s^2]
    const double q_star,          // moisture scale [kg water / kg air]
    const double dtau_du,         // derivative of zonal wind stress w.r.t. the lowest zonal level wind speed [Pa/(m/s)]]
    const double dtau_dv,         // derivative of meridional wind stress w.r.t. the lowest meridional level wind speed [Pa/(m/s)]
    //const double frac_open_sea, // non-seaice fraction - not yet

    In/Out parameters
    double u_flux,                // zonal wind stress [Pa]
    double v_flux,                // meridional wind stress [Pa]

    Output parameters
    double gust,                  // gustiness factor [m/s]

    cosz and solar were already calculated in the F90-side, and thus
    are now input-parameters for Aeolus

    # TODO: Aeolus2 most probably needs "standard" cosz, without radiation weighting
    const double cosz,            // cosine of the zenith angle, per latitude, 1-D array (jcs:jce)
    #const double cosz_radwt,     // cosine of the zenith angle, per latitude, 1-D array (jcs:jce)
    #                             // Note that the SWR code in Aelus needs radiation-weighted cosz!
    const double solar,           // solar irradiation at top of atmosphere, per latitude, 1-D array (jcs:jce)

    To avoid non-linear interference with regular-to-reduced grid
    interpolation, we calculate the upgoing long wave radiation at the
    surface, according to Stefan-Boltzmann-Law, from t_surf already in
    the Fortran interface, and pass it to here.

    double flux_lwr_sf_up,        // upgoing lwr at surface, calculated as in surface_flux()

    // output parameters again
    double flux_sw,               // net shortwave flux (W/m2) at the surface
    double flux_sw_dir,
    double flux_sw_dif,
    double flux_sw_down_vis_dir,
    double flux_sw_down_vis_dif,
    double flux_sw_down_total_dir,
    double flux_sw_down_total_dif,
    double flux_sw_vis,
    double flux_sw_vis_dir,
    double flux_sw_vis_dif,
    double flux_lw,               // DC: downward longwave flux (W/m2) at the surface. The FMS coupler computes the upgoing lwr flux on its own, and then the net LWR flux at surface
   // output parameters from surf_diff_type fields
    double surf_diff_delta_t,     // the increment in temperature in the lowest atmospheric
                                  // layer (((i+1)-(i-1) if atmos model is leapfrog) (K)
                                  // (defined in  gcm_vert_diff_down as the increment computed up
                                  // to this point in model, including effect of vertical
                                  // diffusive flux at top of lowest model level, presumed
                                  // to be modified to include effects of surface fluxes
	                          // outside of this module, then used to start upward
	                          // tridiagonal sweep,
    double surf_diff_dflux_t,     // derivative of the temperature flux at the top of the lowest
                                  // atmospheric layer with respect to the temperature
                                  // of that layer  (J/(m2 K))
    double surf_diff_delta_q,     // - similarly for the increment in specific humidity
                                  //   (non-dimensional  = Kg/Kg)
    double surf_diff_dflux_q,     // - derivative of the flux of specific humidity
                                  //   at the top of the lowest atmospheric layer with respect to
                                  //   the specific humidity of that layer  (--/(m2 K))
    double surf_diff_dtmass,      // dt/mass, where dt = atmospheric time step (sec)
                                  // mass = mass of lowest atmospheric layer (Kg/m2)
    double surf_diff_delta_u,     //
    double surf_diff_delta_v,     //
    If CO2 is defined as tracer in the coupler:
    double surf_diff_delta_co2,   // - similarly for the increment in CO2
                                  //   (non-dimensional  = Kg/Kg)
    double surf_diff_dflux_co2,   // - derivative of the flux of CO2
                                  //   at the top of the lowest atmospheric layer with respect to
                                  //   the specific humidity of that layer  (--/(m2 K))

    """
    global land_frac

    # unpack the argument tuple. Be careful with the ordering!
    Time_seconds, Time_days, dayoftheyear, \
    land_frac, t_surf, albedo, \
    albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif, \
    rough_mom, \
    u_star, b_star, q_star, dtau_du, dtau_dv, \
    u_flux, v_flux, \
    cosz, solar, \
    flux_lwr_sf_up \
    = x

    # TODO: Masoud+Sullyandro: implement this function by invocation of
    # broken-out pieces from (former) Aeolus2_1.py.
    print("TODO: Masoud+Sullyandro: implement aeolus2_update_down()", flush=True)
    # dummy initialisation of output variables
    # as placeholder for the real implementation
    gust = 0.0
    flux_sw = 0.0
    flux_sw_dir = 0.0
    flux_sw_dif = 0.0
    flux_sw_down_vis_dir = 0.0
    flux_sw_down_vis_dif = 0.0
    flux_sw_down_total_dir = 0.0
    flux_sw_down_total_dif = 0.0
    flux_sw_vis = 0.0
    flux_sw_vis_dir = 0.0
    flux_sw_vis_dif = 0.0
    flux_lw = 0.0
    surf_diff_delta_t  = 0.0
    surf_diff_dflux_t  = 0.0
    surf_diff_delta_q = 0.0
    surf_diff_dflux_q = 0.0
    surf_diff_dtmass   = 0.0
    surf_diff_delta_u  = 0.0
    surf_diff_delta_v  = 0.0
    surf_diff_delta_co2 = 0.0
    surf_diff_dflux_co2 = 0.0

    return u_flux, v_flux, \
        gust, \
        flux_sw, \
        flux_sw_dir, flux_sw_dif, \
        flux_sw_down_vis_dir, flux_sw_down_vis_dif, flux_sw_down_total_dir, flux_sw_down_total_dif, \
        flux_sw_vis, flux_sw_vis_dir, flux_sw_vis_dif, \
        flux_lw, \
        surf_diff_delta_t, surf_diff_dflux_t, \
        surf_diff_dtmass, \
        surf_diff_delta_u, surf_diff_delta_v, \
        surf_diff_delta_co2, surf_diff_dflux_co2

# end of aeolus2_update_down()


def aeolus2_update_up(x):
    """
    const int Time_year, Time_month, Time_seconds, // seconds of day
    const int Time_days // days since start of time axis, not day of year!

    All following parameters are 2-D arrays,
    Defined on the local data domain with index range (ics:ice,jcs:jce)

    Input parameter
    const double land_frac,  // fraction amount of land in a grid box
    // two input parameter fields from Surf_diff_type
    const double surf_diff_delta_t,  // the increment in temperature in the lowest atmospheric
                                     // layer (((i+1)-(i-1) if atmos model is leapfrog) (K)
                                     // (defined in  gcm_vert_diff_down as the increment computed up
                                     // to this point in model, including effect of vertical
                                     // diffusive flux at top of lowest model level, presumed
                                     // to be modified to include effects of surface fluxes
	                             // outside of this module, then used to start upward
	                             // tridiagonal sweep,
    const double surf_diff_delta_q,  // - similarly for the increment in specific humidity
                                     //   (non-dimensional  = Kg/Kg)


    Output parameters
    double lprec,  // mass of liquid precipitation since last time step (Kg/m2/s)
    double fprec,  // mass of frozen precipitation since last time step (Kg/m2/s)
    double gust,   // gustiness factor

    Input parameters again
    # u_star, b_star, q_star are not used in EBM atmosphere_up
    const double u_star, // friction velocity
    const double b_star, // bouyancy scale
    const double q_star, // moisture scale
    const int fortytwo   // minimal interface check
    """
    global land_frac

    # unpack the argument tuple. Be careful with the ordering!
    Time_year, Time_month, Time_seconds, \
    Time_days, \
    land_frac, \
    surf_diff_delta_t, surf_diff_delta_q, \
    u_star, b_star, q_star \
    = x

    # TODO: Masoud+Sullyandro: implement this function by invocation of
    # broken-out pieces from (former) Aeolus2_1.py.
    print("TODO: Masoud+Sullyandro: implement aeolus2_update_up()", flush=True)

    # dummy initialisation of output variables
    # as placeholder for the real implementation
    lprec = 0.0
    fprec = 0.0
    gust  = 0.0

    return lprec, fprec, gust

# end of aeolus2_update_up()


# functions for global min/max/sum/mean,
# _with_ and without area weighting.

# Use just global maximum or minimum function
def global_amin(_comm, array):
    localmin = np.amin(array)
    if (_comm is None or _comm.size == 1): return localmin
    return comm.allreduce(localmin, MPI.MIN)

def global_amax(_comm, array):
    localmax = np.amax(array)
    if (_comm is None or _comm.size == 1): return localmax
    return _comm.allreduce(localmax, MPI.MAX)

def global_sum(_comm, array, weights=1.0):
    localsum = (array*weights).sum()
    if (_comm is None or _comm.size == 1): return localsum
    return _comm.allreduce(localsum, MPI.SUM)

def global_mean(_comm, array):
    """
    Unweighted mean over all MPI tasks.
    """
    if (_comm is None or _comm.size == 1): return array.mean()
    localsize = array.size
    globalsize = _comm.allreduce(localsize, MPI.SUM)
    return global_sum(_comm, array)/globalsize

def global_average(_comm, array, weights=None, weightsum=None):
    """
    Weighted mean over all MPI tasks.
    weights=None means weighting by 1.0, as in numpy.average().
    If specified, weightsum is the global sum of weights.
    """
    if (_comm is None or _comm.size == 1): return np.average(array, weights)
    if not weightsum is None:
        return global_sum(_comm, array*weights)/weightsum
    else:
        localsize = array.size
        globalsize = _comm.allreduce(localsize, MPI.SUM)
        return global_sum(_comm, array)/globalsize

def global_count_nonzero(_comm, array):
    localcount = np.count_nonzero(array)
    if (_comm is None or _comm.size == 1): return localcount
    return _comm.allreduce(localcount, MPI.SUM)
