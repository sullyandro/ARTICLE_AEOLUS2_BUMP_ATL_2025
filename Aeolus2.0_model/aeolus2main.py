"""
####################################################################################################
 
  This is the Aeolus2.0 main script.
 
  Author :     Masoud Rostami (rostami@pik-potsdam.de)
  Date:        10.10.2023
  Type:        Dynamical Core
  Institution: Potsdam Institute for Climate Impact Research (PIK)
  Webpage:     www.pik-potsdam.de/en/institute/departments/earth-system-analysis/models/aeolus-2.0
 
 THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 Model description:

    The dynamical core of this atmospheric model is based on the pseudo-spectral moist-convective 
    Thermal Rotating Shallow Water (mcTRSW) model covering the full sphere. The numerical methods 
    used in the model are according to the Dedalus algorithm, which employs spin-weighted 
    spherical harmonics.
    
    Detailed explanations of the model's dynamical core can be found in the following articles. 
    There may be some differences in the implementation of the codes compared to the explanations 
    provided in the articles, which are clarified by the inline comments.
    
    1) Rostami, M., Zhao, B., Petri, S., 2022. On the genesis and dynamics of Madden-Julian 
       oscillation-like structure formed by equatorial adjustment of localized heating. 
       Q. J. R. Meteorol. Soc., 148(749), 3788– 3813, https://doi.org/10.1002/qj.4388.
       
    2) Rostami, M., Severino, L., Petri, S., & Hariri, S., 2023. Dynamics of localized extreme 
       heatwaves in the mid-latitude atmosphere: A conceptual examination. 
       Atmos. Sci. Lett., e1188, https://doi.org/10.1002/asl.1188.
    
####################################################################################################
"""

print()
print('########################################################################################################')
print('#                                                                                                      #')
print('# Starting the Dynamical model core                                                                    #')
print('#                                                                                                      #')
print('# Model name:  Aeolus2.0                                                                               #')
print('# Institution: PIK, Potsdam, Germany                                                                   #')
print('# Webpage:     www.pik-potsdam.de/en/institute/departments/earth-system-analysis/models/aeolus-2.0     #')
print('#                                                                                                      #')
print('########################################################################################################')
print()


import sys
import os
import time ; start_time = time.time()
import json 
import argparse
import math                     as ma
import xarray                   as xr
import numpy                    as np
import scipy.integrate          as integrate #Gives access to the ODE integration package
import scipy.special            as ss
import scipy.io                 as sio
import sphere_wrapper           as sph
import equations_aeolus2main    as eq
import dedalus.public           as de
import warnings
from   datetime                 import timedelta
from   mpi4py                   import MPI
from   dedalus.tools.config     import config # Load config options
from   NetCDFOutput             import NetCDFOutput
import timesteppers

# Suppress FutureWarning messages
warnings.simplefilter(action='ignore', category=FutureWarning)


# --------------------- *** insolation and orbital parameters over the last 5 million years *** -----------
# Note: Consider organizing this section as functions or methods for better readability.
import climlab
from   climlab                  import constants as const
# Importing only necessary functions/classes from climlab
from   climlab.solar.insolation import daily_insolation
from   climlab.solar.orbital    import OrbitalTable
from   climlab.utils.thermo     import pseudoadiabat
#----------------------------------------

# Import user settings from config.py file
from   config import *            # import constants and parameters from the config file
import config as config_ael
varconf = vars(config_ael)
# --------------------- *** End of imports *** ----------------------------------------------------------
#    --------------------- *** User setting & parameterization *** -----------

print()
print('Running --> User setting & parameterization')
print()

'''
varconf['smooth_run']         = smooth_run
varconf['fast_run']           = fast_run
varconf['super_fast']         = super_fast
varconf['restart_run']        = restart_run
varconf['moist_convection']   = moist_convection
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
varconf['n_iterations']       = n_iterations
varconf['delta1']             = delta1
varconf['H0_dim']             = H0_dim     
varconf['g_real']             = g_real    
varconf['U_scale']            = U_scale 
varconf['u_asr']              = u_asr  
varconf['p_level1']           = p_level1  
varconf['p_level2']           = p_level2  
varconf['a']                  = a       
'''
 
print('Nondimensional radius of Earth', a) 
print('Velocity Scale', U_scale)

g0                  = 1.0
second_order_q_eff  = 1                         # 1: includes effect of q_i/h_i on the rhs of h_i and b_i equations.
vaporization_cloud  = 1                         # 1: activates vaporization + precipitable water as clouds for the lower and upper troposphere
uniform_DD          = 0                         # 1: uniform downdraft when downdraft is 1 (ON), 0: nonuniform downdraft
downdraft           = 1                         # 1: active downdraft, 0: no downdraft. With real bottom topo downdraft should be 1

Om                  = a*1.0/2.                  # Omega. It affects parametrization of Coriolis; f=2*Omega*Cos(theta) and scaling 
#sigma_Boltzmann    = 5.671**(-8.)              # the Stefan-Boltzmann constant
L_v                 = 2.5*10.0**(6.0)           # Latent heat of vaporization in J/kg
T0                  = 273.15                    # zero degree in Kelvin
#Co_Cl              = 1.25                      # 1.25, amplification factor of the land–ocean warming contrast, 
                                                # bigger values represent bigger contrast; indeed, Sensible Heat Flux can be merged to this parameter too.
                                                # in improved versions the contrast should be larger for drier land regions and a function of lat.
R_dry               = 287.05                    # R_dry​ is the specific gas constant for dry air [J/(kg.K)].
R_moist             = 461.50                    # R_moist​ is the specific gas constant for water vapor [J\(kg.K)].
tau_bu              = 10.0                      # numerical relaxation time for smooth buoyancy forcing
                                  
# --------------------- *** End of user setting *** -----------------------------------

# Some of the user settings can be overridden via the command line
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--grid',               type=str, help='grid resolution [smooth, fast, super_fast]')
parser.add_argument('-i', '--n_iterations',       type=int, help='total iterations')
parser.add_argument('-d', '--output_folder',      type=str, help='output folder')
parser.add_argument('-o', '--n_output',           type=int, help='data output cadence')
parser.add_argument('--moist_convection',         type=int, help='moist_convection = 0 or 1')
parser.add_argument('--external_forcing',                   help='external_forcing = 0, baroclinic, barotropic')
parser.add_argument('--external_forcing_epsilon',           help='external_forcing_epsilon = 0.1 or 0.2 or 0.3 ...')
args = parser.parse_args()

print(args)

if args.output_folder:
    output_folder            = args.output_folder
    varconf['output_folder'] = output_folder
print('output_folder ', output_folder)    

if args.moist_convection:
    moist_convection            = int(args.moist_convection)
    varconf['moist_convection'] = moist_convection
print('moist_convection ', moist_convection)    

if args.external_forcing:
    external_forcing            = args.external_forcing
    varconf['external_forcing'] = external_forcing
print('external_forcing ', external_forcing)    

if args.external_forcing_epsilon:
    external_forcing_epsilon            = float(args.external_forcing_epsilon)
    varconf['external_forcing_epsilon'] = external_forcing_epsilon
print('external_forcing_epsilon ', external_forcing_epsilon)    


if args.grid:
    if args.grid   == 'smooth':
        smooth_run =  1; fast_run = 0; super_fast = 0;
    elif args.grid == 'fast':
        smooth_run =  0; fast_run = 1; super_fast = 0;
    elif args.grid == 'super_fast':
        smooth_run =  0; fast_run = 0; super_fast = 1;
    else:
        print('Error: grid must be one of [smooth, fast, super_fast] but is ', args.grid)
        sys.exit(42)

print()
print('Running --> MPI configuration')
print() 

print('pid                 ', os.getpid()) 
print('commandline         ', sys.executable, sys.argv) 

# Find MPI rank
comm      = MPI.COMM_WORLD
rank      = comm.rank
size      = comm.size
print('MPI rank            ', rank, ' of ', size)

def print0(arg1='', arg2='', arg3='', arg4='', arg5='', arg6='', arg7='', arg8='', arg9='', arg10=''):
    if rank == 0:
        print(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)

# Setup outputs
if 0 == rank:
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
# make sure that no task attempts to create files inside output_folder
# before it is actually created
comm.Barrier()

# print some information about the execution environment

print0('python version      ', str(sys.version).replace('\n',''))
print0('python installed in ', sys.executable)
print0('running on          ', sys.platform)
print0()
print0('python path is:   \n', str(sys.path).replace('[','').replace(']','').replace(',','\n'))
print0()
print0('float info:       \n', sys.float_info)
print0()
print0('Numpy version       ', np.version.full_version)
print0('Numpy default floating point error handling', np.geterr())
# Sigh. Dedalus provides version info only in its main module,
# but that must not be imported here.
# dedalus.__version__
# dedalus.__path__
print0()
print0('Dedalus public module info:\n', de)
    
# np.seterr(all='raise')

print0()
print0('Running --> Model configurations')
print0()
    
# Define functions for common operations (e.g., global_min, global_max, etc.)
def distance(lon1, lat1, lon2, lat2):
    # return central circle
    dlat   = lat2 - lat1
    dlon   = lon2 - lon1
    q      = np.sin(dlat/2.0)**2.0 + (np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2.0)**2.0))
    output = 2.0 * ma.atan2(np.sqrt(q), np.sqrt(1.0-q))
    return output
# Use just global maximum or minimum function, replace 'a' with your desired variable
def global_amin(comm, a):
    localmin = np.amin(a)
    if (comm == None or comm.size == 1): return localmin
    return comm.allreduce(localmin, MPI.MIN)

def global_amax(comm, a):
    localmax = np.amax(a)
    if (comm == None or comm.size == 1): return localmax
    return comm.allreduce(localmax, MPI.MAX)

def global_sum(comm, a):
    localsum = a.sum()
    if (comm == None or comm.size == 1): return localsum
    return comm.allreduce(localsum, MPI.SUM)

def global_mean(comm, a):
    if (comm == None or comm.size == 1): return a.mean()
    localsize = a.size
    globalsize = comm.allreduce(localsize, MPI.SUM)
    return global_sum(comm, a)/globalsize

def global_count_nonzero(comm, a):
    localcount = np.count_nonzero(a)
    if (comm == None or comm.size == 1): return localcount
    return comm.allreduce(localcount, MPI.SUM)

def save_data(comm, fname, **kwargs):
    datadict = {}
    for name,value in kwargs.items():
        globval = comm.gather(value, root=0)
        if rank == 0:
            datadict[name] = np.hstack(globval)
         
    if rank == 0:
        np.savez(fname, **datadict)
# end of save_data

# these few lines (PERMC_SPEC = None) remove error of timesteppers and is repeated in that script too.
'''
STORE_LU            = config['linear algebra'].getboolean('store_LU') 
PERMC_SPEC          = None #config['linear algebra']['permc_spec']
USE_UMFPACK         = config['linear algebra'].getboolean('use_umfpack')
'''

period              = 2.0*np.pi/Om
gr                  = g0*(B2-1.0)

'''
equivalent_depth = 0
if equivalent_depth==1:
    beta            = 2.0*Om/a
    Ld              = a/5.0
    Lsize           = Ld
    Hr              = 0.2  # H1/H2 lower/upper
    He              = (Ld**4*beta**2.0)/gr
    print('He is : ',He)
    cg1             = np.sqrt(g0*He)
    H0              = He*(1.0+Hr)**2.0/Hr
    H1              = H0*Hr/(1.0+Hr)
    H2              = H0/(1.0+Hr)
    spara           = 4
else:# H0=1 is total troposphere/atmosphere
'''
#if topography == 1:
#   delta1          = 0.60#0.35#0.57 # H1/H0
#else:
#   delta1          = 0.5#0.35
delta2              = 1.0-delta1
H0                  = 1.0
H1                  = delta1*H0
H2                  = delta2*H0
cg1                 = np.sqrt(g0*H0)

# Discretization parameters
L_dealias           = 3/2
if smooth_run == 1:
    L_max           = 255       # spherical harmonic order
    S_max           = 3         # spin order (leave fixed)
    step            = 1.        # minutes, numerical time step
    Nxx             = 768.
    Nyy             = 384.    
    theta_basis = de.Fourier('theta', L_max+1, interval=(0,np.pi), dealias=L_dealias)
if fast_run   == 1:
    L_max           = 127       # spherical harmonic order
    S_max           = 3         # spin order (leave fixed)
    step            = 2.5       # minutes, numerical time step, the model is convergent at step=5
    Nxx             = 384.
    Nyy             = 192.    
    theta_basis = de.Fourier('theta', L_max+1, interval=(0,np.pi), dealias=L_dealias)   
if super_fast == 1:
    L_max           = 63        # spherical harmonic order
    S_max           = 3         # spin order (leave fixed) 
    step            = 5.0       # 5 is ok., tested, numerical time step
    Nxx             = 192.
    Nyy             = 96. 
    theta_basis     = de.Fourier('theta', L_max+1, interval=(0,np.pi), dealias=L_dealias)  


# Integration parameters
dt                  = step*60.*period/86400.  # Equal to one minute numerical time step
                             
if args.n_iterations:
    n_iterations = args.n_iterations
#n_output     = np.around(a*58, 0)# np.around(2*a*58, 0) data output cadence
n_output     = np.around(((1./T_scale)/dt)/cadence, 0)# np.around(2*a*58, 0) data output cadence

if rank == 0: print('Number of iterations per output = ', ((1./T_scale)/dt)/cadence)

if args.n_output:
    n_output        = args.n_output


print0('dt                  ', dt)
print0('n_iterations        ', n_iterations)
print0('n_output            ', n_output)


if topography == 1:
    t_p_inv         = 1./(35.0*(dt/step))   # 1./(30.*(dt/step))  # convective adjustment/relaxation time scale, 1/(condensation relaxation time) 
                                            # bigger than 30, leads to higher condensation on the equator
    gamma_NC1       = 5.0                   # 3.5 coefficient of the Newtonian Cooling force (can be different from gamma and for each layer)
    gamma_NC2       = 5.0                   # 3.5
    gamma_RT1       = 1.0                   # 1.0001, coefficient of the Radiative Transfer flux in the lower layer (can differ from gamma_* and vary for each layer)
    gamma_RT2       = 1.0                   # 1.0004, coefficient of the Radiative Transfer flux in the upper layer (can differ from gamma_* and vary for each layer)
                                            # *** Notion: Values should be lower than 1 and not bigger than 1.004 that leads to an abrupt variation of h1 and h2 and
                                            # unrealistic results. Lower values in Jan. --> bigger u1 and u2
    gamma_topo      = 8.0                   # 6.5, bigger = reducing the bottom topographic effects  
    nu              = step*0.00009          # 0.00009 Numerical Laplacian damping for momentum and buoyancy equations to preclude small-scale convective instabilities
    kappa           = 1.1*nu                # the same term for h  
else:
    t_p_inv         = 1./(30.*(dt/step))    # 1/(condensation relaxation time)  
    gamma_NC1       = 1.0 #2.5              # coefficient of the Newtonian Cooling force (can be different from gamma and for each layer)
    gamma_NC2       = 1.0 #3.5
    gamma_RT1       = 1.0                   # coefficient of the Radiative Transfer flux in the lower layer (can differ from gamma_* and vary for each layer)
    gamma_RT2       = 1.0                   # coefficient of the Radiative Transfer flux in the upper layer (can differ from gamma_* and vary for each layer)
                                            # gamma_RT should not be much bigger than 1. 
    gamma_topo      = 1.0                   # 1.7# bigger = reducing the bottom topographic effects  
    nu              = step*0.0001           # 2.e-5 Newtonian viscosity/diffusivity/thermoconductivity for momentum and buoyancy equations to preclude small-scale
                                            #  convective instabilities
    kappa           = 0.7*nu                # the same term for h    
    coeff_RT1        = 1.0
    coeff_RT2        = 1.0
    
t_v_inv             = (1.0/20.0)*t_p_inv    # 1/(vaporization relaxation time)
t_drag_inv          = 0.                    # 1/tau_u, just for TQG adjustment 

if vaporization_cloud == 1:
    ivap            = 1.0                   # Flag indicating vaporization
else:
    ivap            = 0.0

           
if moist_convection == 1:
    i_mc            = 1.0
    if topography == 1:                     # Flag for enabling topography effects  
        Qs1         = 0.0145#0.9*H1 # 0.01         # physically it can be fixed to 0.9 and H_i be included in ep_i but the results will be the same
        Qs2         = 0.9*H2                # 0.01 # 0.04  
        Q01         = 0.                    # Qs1*(1.0-0.015)
        Q02         = 0.                    # Qs2*(1.0-0.015)
        max_q       = 2.0*0.02#Qs1              # 1.03*Qs1 # 5.0*(Qs1-Q01)
        gamma       = 0.70                  # 0.45# Lower = more intense condensation with higher latent heat release
        ep1         = 0.1*(1.0 - gamma)     # 2, Condensation effciency parameter, e.g. =Lv/Cp/d_theta, nondimensional order of CC1*ep1 is about 1e-2
        t_r_b1      = 1./26.0*dt#1./10.0*dt#0.*1.0/(22.0*T_scale) #0.#  Relaxation to background SST or b10
        t_r_b2      = 1./26.0*dt#1./10.0*dt#0.*1.0/(27.0*T_scale) # 0.# Bigger relaxation leads to strengthening the insolation. 
        wcritical   = 0.0002                # 0.0002 threshold value of precipitation.      
    else:
        Qs1         = 0.9*H1 # 0.04         # physically it can be fixed to 0.9 and H_i be included in ep_i but the results will be the same
        Qs2         = 0.9*H2 # 0.04   
        Q01         = Qs1*(1-0.01)          # Close to saturation value
        Q02         = Qs2*(1-0.01)
        t_e_inv     = 15. # 3.0 # 4.0       # for conceptual cases. Indeed, it can be time independent, but also a function of SST. 0.7 is a stable parameter
        coeff_v     = 0.1*t_e_inv           # 0.1 efficiency/effect of the surface velocity upon sea surface evaporation        
        gamma       = 0.1 # 0.3
        t_r_b1      = 0.*1.0/(22.0*T_scale) # 1./10*dt#0.#      /4*60*dt # in literature is more than 25 days for SST
        t_r_b2      = 0.*1.0/(22.0*T_scale) # 1./10*dt#0.#     /4*60*dt
        ep1         = 1.0 - gamma           # condensation efficiency =Lv/Cp/d_theta
        wcritical   = 0.01#0.0001           # threshold value of precipitation.      
    upward          = 0.3                   # upward draft of humidity from the lower troposphere to the upper layer
    t_r_u1          = 0.                    # applies just for getting TQG adjustment, otherwise is zero
    t_r_u2          = t_r_u1
elif unparallel_TQG_adj == 1:
    i_mc            = 0.
    iradN           = 0.
    Qs1             = 0.9*H1
    Qs2             = 0.9*H2  
    max_q           = 1.06*Qs1
    gamma           = 0.1 # 0.3 
    ep1             = 1.0*(1.0 - gamma)     # condensation efficiency =Lv/Cp/d_theta
    Q01             = Qs1*(1.0-0.01)
    Q02             = Qs2*(1.0-0.01)   
    t_e_inv         = 0.0
    
    #--------------- *** freezing u2,b1,b2 (just three of them should be freezed to get TQG balanced state) *** --------------
                                            # The best choice is freezing u2,b1,b2, if u2 is an active layer.
    t_r_b1          = 0.*1./(4*dt)          # freezing b1, applies just for getting TQG adjustment, otherwise is zero
    t_r_b2          = 0.*1./(4*dt)          # freezing b2, applies just for getting TQG adjustment, otherwise is zero
    t_r_u1          = 0.*1./(4*dt)          # freezing u1, applies just for getting TQG adjustment, otherwise is zero
    t_r_u2          = 0.*1./(4*dt)          # freezing u2, applies just for getting TQG adjustment, otherwise is zero
    #------------------------------------------ end of freezing ------------------------------

    upward          = 0.                    # upward draft of humidity from the lower troposphere to the upper layer
    wcritical       = 0.005                 # threshold value of precipitation.  
else: # dry case
    i_mc            = 0.
    Qs1             = 0.9*H1
    Qs2             = 0.9*H2   
    Q01             = Qs1*(1-0.01)
    Q02             = Qs2*(1-0.3)   
    t_e_inv         = 0.0
    t_r_b1          = 0. # 1./(5*dt)        # 1/4*60*dt#, restoring b1
    t_r_b2          = 0. # 1./(5*dt)        # 1/4*60*dt#, restoring b2
    t_r_u1          = 0. # 1./(4*dt)        # applies just for getting TQG adjustment, otherwise is zero
    t_r_u2          = t_r_u1
    ep1             = 0. 
    upward          = 0.                    # upward draft of humidity from the lower troposphere to the upper layer
    wcritical       = 0.005                 # threshold value of precipitation. 
   
if topography == 1:
    itopo           = 1.0
else:
    itopo           = 0.0
   
#if second_order_q_eff == 1:
#    iq_o_h          = 0.01 # 20.744 # 1.0/200. # re-scaling for arbitrary value of: q/h[dimensional]=[L_v.g/(C_p.theta_s)].q*/h*, where * indicates nondimensional value. q_real/q_model=45
#else:
#    iq_o_h          = 0.0 

if Newtonian_cooling == 1:
    iradN           = 1.0
    iradh           = 0.0
    t_r_inv         = 1.0/(90.0*T_scale)    # 1.0/(50.0*T_scale)# Thermal radiation relaxation time/ Newtonian cooling relaxation time ~ 90 days, if it relaxes to initial condition, otherwise it can be bigger.    
else:
    iradN           = 0.0   
    iradh           = 0.0 # 1.0: for choosing another cooling option
    t_r_inv         = 0.0   
if Radiative_Transfer == 1:
    iRT             = 1.0
else:
    iRT             = 0.0   
   
if downdraft == 1:
    iDD             = 1.0
else:
    iDD             = 0.0     

# print('theta_basis', theta_basis) 
# Make domain
lamda_basis         = de.Fourier('lamda', 2*(L_max+1), interval=(0,2*np.pi), dealias=L_dealias)
domain              = de.Domain([lamda_basis,theta_basis], grid_dtype=np.float64, comm=comm)

# set up sphere
m_start             = domain.distributor.coeff_layout.start(1)[0]
m_len               = domain.distributor.coeff_layout.local_shape(1)[0]
m_end               = m_start + m_len - 1
N_theta             = int((L_max+1)*L_dealias)


print0()
print0('L_max               ', L_max)
print0('S_max               ', S_max)
print0('L_dealias           ', L_dealias)
print0('m_start             ', m_start)
print0('m_len               ', m_len)
print0('m_end               ', m_end)
print0('N_theta             ', N_theta)

# def __init__(self,L_max,S_max=0,N_theta=None,m_min=None,m_max=None): The maximum spherical harmonic degree is L max , and the azimuthal order m ranges 0 ≤ m ≤ L max .

S                   = sph.Sphere(L_max,S_max,N_theta=N_theta,m_min=m_start,m_max=m_end)
#print("S grid, weights, sin_grid, cos_grid", S, S.grid, S.weights, S.sin_grid, S.cos_grid) 
lamda               = domain.grids(L_dealias)[0]
theta_slice         = domain.distributor.grid_layout.slices(domain.dealias)[1]
theta_len           = domain.local_grid_shape(domain.dealias)[1]
theta_global        = S.grid
theta               = S.grid[theta_slice].reshape((1,theta_len))
#phi                = np.pi/2.0-theta
phi                 = theta-np.pi/2.0

# print('theta_len        ', theta_len)
# print('theta_global.size', theta_global.size)
# print('theta_slice.start', theta_slice.start)
# print('theta_slice.stop ', theta_slice.stop)
# print('theta_slice.step ', theta_slice.step)
# print('theta_slice      ', theta_slice)
# print('theta global     ', theta_global)
# print('theta            ', theta) 
# print('theta in deg     ', theta/(2*np.pi)*360.0) 
# print('phi              ', phi)

# Sigh. The current input topographies are organized as -180 to 180 deg longitude, but lamda runs 0 to 2pi. Thus we need to rotate by -180deg to obtain a correct axis for NetCDF. 
# Sigh. This results in Y running N to S, while FMS wants coordinates running S to N. This is consistent with the current input topographies. ncview and ferret both revert flip the latitudes automagically, with more or less warnings printed out. 

#lons               = lamda/(2.*np.pi)*360.0
lons                = lamda/(2.*np.pi)*360.0 - 180.0

lats                = phi/(2.*np.pi)*360.0
#lats               = theta/(2*np.pi)*360.0 - 90.0



# Indices of the local domain boundaries in the global domain. Naming convention taken from GFDLs FMS coupler. i/j/k denotes lon/lat/vertical direction, c/d for compute resp. data domain (i.e. without/with halos), s/e for start/end 
# sigh. "is" is a reserved word in python
ics                 = 0                 
ice                 = lons.size
jcs                 = theta_slice.start
jce                 = theta_slice.stop

#print('Lon size, Lat size: ', lons.size, lats.size)
#print('lons',    lons)
#print('lats',    lats)
#print('Lons shape, Lats shape: ', lons.shape, lats.shape) 

gridinfo = {
    'nlons_global' : lons.size,
    'nlats_global' : theta_global.size,
    'lons'         : lons[:,0],
    'lats'         : lats[0,:],
    'ics'          : 0,
    'ice'          : lons.size,         # according to python conventions, one after end
    'jcs'          : theta_slice.start,
    'jce'          : theta_slice.stop,  # according to python conventions, one after end
}

ncout = NetCDFOutput(comm, gridinfo, os.path.join(output_folder, 'Aeolus2-output.nc'), comment='Testing NetCDFOutput.py')

# Here we must use the names of the axes as defined in the NetCDF file
axes_ty  = (ncout.timeaxisname, ncout.lataxisname,)                    # time series of zonal-mean values
axes_yx  = (ncout.lataxisname,  ncout.lonaxisname,)                    # static maps 
axes_tyx = (ncout.timeaxisname, ncout.lataxisname, ncout.lonaxisname,) # time series of maps

#--------------- idealized temperature profile -------
def generate_idealized_temp_profile(T_s, plevs, Tstrat=200):
    '''
    Generation of an idealized temperature profile with specified T_s and Tstrat
    '''
    solution                    = integrate.odeint(pseudoadiabat, T_s, np.flip(plevs))
    temp                        = solution.reshape(-1)
    temp[np.where(temp<Tstrat)] = Tstrat
    return np.flip(temp)                   # need to be reinverted with respect to the pressure axis
# end of generate_idealized_temp_profile() 

    
def make_idealized_column(T_s, num_lev=100, Tstrat=200):
    state            = climlab.column_state(num_lev=num_lev, num_lat=1)   # Set up a column state
    plevs            = state['Tatm'].domain.axes['lev'].points            # Extract the pressure levels
    state['Ts'][:]   = T_s                                      
    state['Tatm'][:] = generate_idealized_temp_profile(T_s=T_s, plevs=plevs, Tstrat=Tstrat)
    return state   
# end of make_idealized_column() 

#------------- *** Begining of insolation *** --------    
# ----- *** insolation for millenial scale, present, or past *** ------ : 
orb_0               = OrbitalTable.interp(kyear=0)                 # present-day orbital parameters
orb_10              = OrbitalTable.interp(kyear=-10)               # orbital parameters for 10 kyrs before present
days                = np.linspace(0, const.days_per_year, 365)
insol               = daily_insolation(lats[0,:], days, orb_0)     
#insol              = daily_insolation(-lats[0,:], days, orb_0)             
if Radiative_Transfer==1:
    pass
else:
    insol           = 0.*insol  

insol_annual_mean = np.sum( np.mean( insol, axis=1 ) * np.cos( np.deg2rad(lats[0,:]) ) ) / np.sum( np.cos( np.deg2rad( lats[0,:])))
print0()
print0('insolation mean     ', insol_annual_mean, '(the area-weighted global, annual average of insolation)')
print0('insol shape         ', insol.shape)
# print('insol [day1]:', insol[:,175])    
            
#--------- millenial scale with orbital params ---------
'''
    insol           = daily_insolation( lats[0,:], days, orb) # for having millenial scale insolation 
    Qann            = np.mean(insol, axis=1)  # time average over the year
    #print('Qann.shape', Qann.shape)
    #print('Lat:', lats[0,:])
    Qglobal         = np.empty_like( kyears )
    for n in range( kyears.size ):   # global area-weighted average
        Qglobal[n]  = np.sum( Qann[:,n] * np.cos( np.deg2rad(lats[0,:]) ) ) / np.sum( np.cos( np.deg2rad(lats[0,:])))
        #print(Qglobal.shape)    
    #print('rows and columns of lats',len(lats), len(lats[0]))
    #print('rows and columns of insol',len(insol), len(insol[0])) 
    print('size of insol',insol.size)
    # To set up the model with different orbital parameters: orb = {'ecc':0., 'obliquity':90., 'long_peri':0.} 
    orb_0           = OrbitalTable.interp(kyear=0)     # present-day orbital parameters
    orb_10          = OrbitalTable.interp(kyear=-10)   # orbital parameters for 10 kyrs before present
    orb_23          = OrbitalTable.interp(kyear=-23)   # 23 kyrs before present
    Q_0             = daily_insolation( lats[0,:], days, orb_0 )
    Q_10            = daily_insolation( lats[0,:], days, orb_10 )   # insolation arrays for each of the three sets of orbital parameters
    Q_23            = daily_insolation( lats[0,:], days, orb_23 )
    Qdiff           = Q_10 - Q_23
    print('The area-weighted global average of the difference:',np.average(np.mean(Qdiff,axis=1), weights=np.cos(np.deg2rad(lats[0,:]))))
'''   
#------------- *** End of insolation *** --------    


print0()
print0('Running --> sph.TensorField for variables')
print0()

# def __init__(self, rank, S, domain):
# rank 0: scalar fields (rank 0), e.g., density, temperature, pressure, the divergence of the velocity
# rank 1: vector fields (rank 1), e.g., velocities, magnetic fields, temperature gradient
# rank 2: tensor fields         , e.g., the strain rate, Maxwell stress; and higher order tensors

u1                  = sph.TensorField(1, S, domain)
h1                  = sph.TensorField(0, S, domain)
u2                  = sph.TensorField(1, S, domain)
h2                  = sph.TensorField(0, S, domain)
q1                  = sph.TensorField(0, S, domain)
b1                  = sph.TensorField(0, S, domain)
b2                  = sph.TensorField(0, S, domain)
u10                 = sph.TensorField(1, S, domain)
u20                 = sph.TensorField(1, S, domain)
b10                 = sph.TensorField(0, S, domain)
b20                 = sph.TensorField(0, S, domain)
h10                 = sph.TensorField(0, S, domain)
q2                  = sph.TensorField(0, S, domain)
w1                  = sph.TensorField(0, S, domain)
w2                  = sph.TensorField(0, S, domain)

Du1                 = sph.TensorField(2, S, domain)
Du2                 = sph.TensorField(2, S, domain)
Db1                 = sph.TensorField(1, S, domain)
Db2                 = sph.TensorField(1, S, domain)
Dh1                 = sph.TensorField(1, S, domain)
Dh10                = sph.TensorField(1, S, domain)
Dh2                 = sph.TensorField(1, S, domain)
uh1                 = sph.TensorField(1, S, domain)
uh2                 = sph.TensorField(1, S, domain)
uq1                 = sph.TensorField(1, S, domain)
uq2                 = sph.TensorField(1, S, domain)
#ub1                = sph.TensorField(1, S, domain)
#uhf1               = sph.TensorField(1, S, domain)
#uqf1               = sph.TensorField(1, S, domain)
#uqf2               = sph.TensorField(1, S, domain)
#ubf1               = sph.TensorField(1, S, domain)
uw1                 = sph.TensorField(1, S, domain)
uw2                 = sph.TensorField(1, S, domain)
#uwf1               = sph.TensorField(1, S, domain)
#uwf2               = sph.TensorField(1, S, domain)

divuh1              = sph.TensorField(0, S, domain)
divuh2              = sph.TensorField(0, S, domain)
divuq1              = sph.TensorField(0, S, domain)
divuq2              = sph.TensorField(0, S, domain)
divuw1              = sph.TensorField(0, S, domain)
divuw2              = sph.TensorField(0, S, domain)
#divub1             = sph.TensorField(0, S, domain)
#divuhf1            = sph.TensorField(0, S, domain)
#divuqf1            = sph.TensorField(0, S, domain)
#divuwf2            = sph.TensorField(0, S, domain)
#divubf1            = sph.TensorField(0, S, domain)


u1_rhs              = sph.TensorField(1, S, domain)
h1_rhs              = sph.TensorField(0, S, domain)
u2_rhs              = sph.TensorField(1, S, domain)
h2_rhs              = sph.TensorField(0, S, domain)
q1_rhs              = sph.TensorField(0, S, domain)
b1_rhs              = sph.TensorField(0, S, domain)
b2_rhs              = sph.TensorField(0, S, domain)
q2_rhs              = sph.TensorField(0, S, domain)
w1_rhs              = sph.TensorField(0, S, domain)
w2_rhs              = sph.TensorField(0, S, domain)

Ev                  = sph.TensorField(0, S, domain) # Evaporation
Oc                  = sph.TensorField(0, S, domain) # Ocean
CC1                 = sph.TensorField(0, S, domain) # Condensation/Latent heat release forcing/ep
DD1                 = sph.TensorField(0, S, domain) # Downdraft+bias correction due to condensation in the lower layer, from upper troposphere to the lower troposphere
DDr                 = sph.TensorField(0, S, domain) # Downdraft without bias correction due to condensation in the lower layer, from upper troposphere to the lower troposphere
VV2                 = sph.TensorField(0, S, domain) # Vaporization in the upper layer
VV1                 = sph.TensorField(0, S, domain) # Vaporization in the lower layer
Prec1               = sph.TensorField(0, S, domain) # Precipitation in the lower layer
#speed              = sph.TensorField(0, S, domain) # velocity magnitude
R2                  = sph.TensorField(1, S, domain)
if Newtonian_cooling == 1:
   Rad1             = sph.TensorField(0, S, domain) # Newtonian cooling for the lower layer
   Rad2             = sph.TensorField(0, S, domain) # Newtonian cooling for the upper layer
RT1                 = sph.TensorField(0, S, domain) # RRTM_LW/SW, lower layer
RT2                 = sph.TensorField(0, S, domain) # RRTM_LW/SW, upper layer
#q_o_h1             = sph.TensorField(0, S, domain) # q1/h1
#q_o_h2             = sph.TensorField(0, S, domain)
#Press1             = sph.TensorField(0, S, domain) # Pressure Lower Layer
#Press2             = sph.TensorField(0, S, domain) # Pressure Upper Layer


state_vector        = eq.StateVector(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)
RHS                 = eq.StateVector(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)
timestepper         = timesteppers.SBDF2(eq.StateVector, u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)


print0()
print0('Running --> Printing coordenate variables')
print0('lons  --> {0:10} [{1:7.3f} \t   {2:7.3f}\t... {3:7.3f} \t  {4:7.3f}]'.format(str( lons.shape),  lons[0,0],  lons[1,0],  lons[-2,0],  lons[-1,0]))
print0('lats  --> {0:10} [{1:7.3f} \t   {2:7.3f}\t... {3:7.3f} \t  {4:7.3f}]'.format(str( lats.shape),  lats[0,0],  lats[0,1],  lats[0,-2],  lats[0,-1]))
print0('lamda --> {0:10} [{1:7.3f} \t   {2:7.3f}\t... {3:7.3f} \t  {4:7.3f}]'.format(str(lamda.shape), lamda[0,0], lamda[1,0], lamda[-2,0], lamda[-1,0]))
print0('phi   --> {0:10} [{1:7.3f} \t   {2:7.3f}\t... {3:7.3f} \t  {4:7.3f}]'.format(str(  phi.shape),   phi[0,0],   phi[0,1],   phi[0,-2],   phi[0,-1]))
# lons  --> (768, 1)   [-180.000 -179.531   ... 179.062      179.531]
# lats  --> (1, 384)   [89.642     89.177   ... -89.177      -89.642]
# lamda --> (768, 1)   [ 0.000      0.008   ...  6.267         6.275]
# phi   --> (1, 384)   [ 1.565      1.556   ... -1.556        -1.565]

#####################################################################
#########################      BUMP     #############################
#####################################################################
print0()
print0('Running --> Bump forcing creation')

# Here you can define an arbitrary forcing

# Bump on Atlantic Ocean   (sh_E = -3.85 / sh_N = 0.36)

# Bump on Sahara Desert    (sh_E = -3.1  / sh_N = 0.36)

# Bump on Turkestan Desert (sh_E = -2.0  / sh_N = 0.78)

if external_forcing != 0:
    
    asr             = 1             # 1.8  # aspect ratio of analytic perturbation: lower--> y elongated, higher--> X-elongated
    if topography == 1:
        sh_E        = -3.85          # shift to east
    else:   
        sh_E        = 0.            # -1.1 #-1.2 # shift to east 
    sh_N            = 0.36          # -0.7 # 0. # positive value = shift to North
    spara           = 4.0           # 4     
    psi             = phi*lamda*0
    r_nd            = phi*lamda*0
    Lsize           = 0.25
    for i in range(lamda.size):
        for j in range(phi.size):
            la, ph  = lamda[i,0], phi[0,j]
            if(la > np.pi):
                la -= 2.0*np.pi
            # x,y = radius*np.cos(ph)*np.minimum(la,2*np.pi-la), radius*ph  r_nd = np.sqrt(x**2+y**2)/Lsize
            # r_nd[i,j] = distance(0,0,la/asr,ph)/Lsize, Lsize low --> compacted
            r_nd[i,j]   = distance(0+sh_E,0+sh_N,la/asr,ph)/Lsize # for zonally elongated bumps       
            # r_nd[i,j] = distance(0+sh_E,0+sh_N,la,ph/asr)/Lsize # for meridionally elongated bumps 
            r_low = r_nd[i,j]**spara/2.
            if r_low > 1.e4:
                r_low   = np.inf
            zeta        = 1./spara+1./2.
            result, err = integrate.quad(lambda t: np.exp(-t) * t**(zeta-1), r_low, zeta)
            psi[i,j]    = result
            if(r_nd[i,j] == 0):
                # print0('psi[i,j]',result )
                # print0('psi[i,j]',psi[i,j])
                pass
    epsilon   = external_forcing_epsilon # from config.py
    # epsilon = 0.21 # 0.4 for h1 , 0.16 for positive b1 anomaly for aquaplanet, 0.35 with real topo, 0.7 for negative h1 anomaly,
    # epsilon = 0.22 is the highest before model explod.
    #h1['g']  = - epsilon * H1* np.sqrt(2.0*np.exp(1.0)) * 2.0**(1./spara)/spara * psi
    #b1['g']  =   epsilon * H1* np.sqrt(2.0*np.exp(1.0)) * 2.0**(1./spara)/spara * psi
    psi_b     = H1* np.sqrt(2.0*np.exp(1.0)) * 2.0**(1./spara)/spara * psi
    f1_ext    = epsilon * (psi_b - global_amin(comm, psi_b))
    #b10['g'] = b1['g']     #0.*h1['g'] 
    #b20['g'] = 0.*b1['g']  #0.*h1['g'] 
    #b10['g'] = 0.*h1['g'] 
    # print('size of lambda and theta:', lamda.size, phi.size)

# b1_bump = epsilon * H1* np.sqrt(2*np.exp(1)) * 2**(1./spara)/spara * psi  # Bump perturbation 
# b1_bump = np.where(b1_bump <= -0.01, 0, b1_bump) 

# b2_bump = epsilon * H2* np.sqrt(2*np.exp(1)) * 2**(1./spara)/spara * psi  # Bump perturbation 
# b2_bump = np.where(b2_bump <= -0.01, 0, b2_bump) 
        
#####################################################################
#####################################################################
#####################################################################

print0()
print0('Running --> Reading topograph and Ocean cover')
print0()

if smooth_run == 1:
    print0('smooth_run reading h10 and Oc from Topo_non_unigrid.mat') 
    h10['g'][0]     = (delta1)*itopo*((sio.loadmat('Topo_non_unigrid.mat')['h_topo_smooth'][ics:ice,jcs:jce]))/(gamma_topo*g_real*H0_dim)
    Oc['g'][0]      =                  sio.loadmat('Topo_non_unigrid.mat')['Ocean'][ics:ice,jcs:jce]
    
if fast_run   == 1: 
    print0('fast_run reading bottom topography (h10) and Ocean (Oc) from Topo_nonunigrid_fast.mat') 
    h10['g'][0]     = (delta1)*itopo*((sio.loadmat('Topo_nonunigrid_fast.mat')['h_topo_fast'][ics:ice,jcs:jce]))/(gamma_topo*g_real*H0_dim)
    Oc['g'][0]      =                  sio.loadmat('Topo_nonunigrid_fast.mat')['Ocean_fast'][ics:ice,jcs:jce]   #size: 384 by 192
    
if super_fast == 1: 
    print0('Super_fast run reading bottom topography (h10) and Ocean (Oc) from Topo_nonunigrid_supfast.mat') 
    h10['g'][0]     = (delta1)*itopo*((sio.loadmat('Topo_nonunigrid_supfast.mat')['h_topo_supfast'][ics:ice,jcs:jce]))/(gamma_topo*g_real*H0_dim)
    Oc['g'][0]      =                  sio.loadmat('Topo_nonunigrid_supfast.mat')['Ocean_supfast'][ics:ice,jcs:jce]   #size: 192 by 96
   
# print('size of lambda      ', lamda.size)
# print('size of theta       ', phi.size)


'''
# introducing analitic forcing for potential temperature or thickness
psib                = phi*lamda*0
for i in range(lamda.size):
    for j in range(phi.size):
        th0         = phi[0,j]
        th2         = (th0+np.pi/2)*180/np.pi
        steep       = 0.007
        y0          = 90
        p1          = 18
        result      = (mp.sech(steep*(th2-y0)))**p1
        psib[i,j]   = 0.0013*result
b1['g']             = psib
b2['g']             = psib
# b10['g']          = psib
# b20['g']          = psib
#------------ end of analytic curve fit of b_i ---------------
'''

file_num            = 1
t                   = 0


# ------ initial conditions for unparallel TQG adjustment -------
if unparallel_TQG_adj == 1:

    #-------- for initial adjustment of desired data ----------------------
    '''
    print('reading initial conditions for unparallel TQG adjustment from Jan1980.mat')   
    b1['g'][0]      =     sio.loadmat('Jan1980.mat')['b1']
    b2['g'][0]      =     sio.loadmat('Jan1980.mat')['b2']
    b10['g'][0]     = 0.0*sio.loadmat('Jan1980.mat')['b1']
    #h10['g'][0]    =    (sio.loadmat('Topo_non_unigrid.mat')['h_topo_smooth']*itopo)/H0_dim
    b20['g'][0]     = 0.0*sio.loadmat('Jan1980.mat')['b1']   
    h1['g'][0]      = 0.0*sio.loadmat('Jan1980.mat')['b1']
    h2['g'][0]      = 0.0*sio.loadmat('Jan1980.mat')['b1']
    u1['g'][1]      =     sio.loadmat('Jan1980.mat')['U1']
    u2['g'][1]      =     sio.loadmat('Jan1980.mat')['U2']
    u1['g'][0]      =     sio.loadmat('Jan1980.mat')['V1']
    u2['g'][0]      =     sio.loadmat('Jan1980.mat')['V2']   
    u10['g'][1]     = 0.0*sio.loadmat('Jan1980.mat')['U1']
    u20['g'][1]     = 0.0*sio.loadmat('Jan1980.mat')['U1']
    '''
    
    # Construct the file name using the month and year variable
    file_name  = f'{month_year}.mat'
    # Print a message indicating the file being read
    print0('Reading initial conditions for unparallel TQG adjustment from {}'.format(file_name))
    # Load data from the specified file
    b1['g'][0]      =       sio.loadmat(file_name)['b1']
    b2['g'][0]      =       sio.loadmat(file_name)['b2']
    b10['g'][0]     = 1.0 * sio.loadmat(file_name)['b1']
    b20['g'][0]     = 1.0 * sio.loadmat(file_name)['b2']
    h1['g'][0]      = 0.0 * sio.loadmat(file_name)['b1']
    h2['g'][0]      = 0.0 * sio.loadmat(file_name)['b1']
    u1['g'][1]      =       sio.loadmat(file_name)['U1']
    u2['g'][1]      =       sio.loadmat(file_name)['U2']
    u1['g'][0]      =       sio.loadmat(file_name)['V1']
    u2['g'][0]      =       sio.loadmat(file_name)['V2']
    u10['g'][1]     = 1.0 * sio.loadmat(file_name)['U1']
    u20['g'][1]     = 1.0 * sio.loadmat(file_name)['U1']



print0('max & min h1        ', global_amax(comm, h1['g']), global_amin(comm, h1['g']))
print0('H1, H2, & H0        ', H1, H2, H0)
print0('min moist enthalpy  ', global_amin(comm, H1-ep1*Q01-ep1*q1['g']), '(layer 1 at time = 0)')


# initialization or restart runs through reading of an old output and continue from then:
if restart_run == 1:
   
    # old_num  --> this name/number should be consistent with your existing file
    if os.path.exists('{}/{}.npz'.format(output_folder, start_file)):
        old_num     = int(start_file.split('_')[-1].replace('.npz','')) # getting number from output/output_"0000N".npz
    else:
        old_num     = 0                                 # if starting from Obs data output_1.npz, then old_num starts in 0
    file_num        = old_num+1
    restart_name    = '{}/{}.npz'.format(input_folder, start_file)
    print0('reading restart file ', restart_name)
    file            = np.load(restart_name)    
    print0('h1 shape ', file['h1'].shape)
    b_eff           = 1.0  # background_strength
    h1['g'][0]      =    file['h1'][ics:ice,jcs:jce]
    h2['g'][0]      =    file['h2'][ics:ice,jcs:jce]
    mean_h1         = 0. #h1['g'][0].mean()
    mean_h2         = 0. #h2['g'][0].mean()
    h1['g'][0]      =  ((file['h1']-mean_h1)*b_eff)[ics:ice,jcs:jce]    # - (epsilon * H1* np.sqrt(2*np.exp(1)) * 2**(1./spara)/spara * psi)    # adding perturbation
    h2['g'][0]      =  ((file['h2']-mean_h2)*b_eff)[ics:ice,jcs:jce]
    u1['g'][0]      =    file['u1th'][ics:ice,jcs:jce]
    u1['g'][1]      =   (file['u1ph']*b_eff)[ics:ice,jcs:jce]
    u2['g'][0]      =    file['u2th'][ics:ice,jcs:jce]
    u2['g'][1]      =   (file['u2ph']*b_eff)[ics:ice,jcs:jce]
    q1['g'][0]      =    file['q1'][ics:ice,jcs:jce]                    # Q01
    b1['g'][0]      =    file['b1'][ics:ice,jcs:jce]                    # + (epsilon * H1* np.sqrt(2*np.exp(1)) * 2**(1./spara)/spara * psi)    # adding perturbation
    # b1['g'][0]    = b1_bump                                           # adding perturbation
    # b2['g'][0]    = b2_bump                                           # adding perturbation
    b2['g'][0]      =    file['b2'][ics:ice,jcs:jce]
    q2['g'][0]      =    file['q2'][ics:ice,jcs:jce]
    b10['g'][0]     = 1.*file['b1'][ics:ice,jcs:jce]                    # b_eff*sio.loadmat('input_b_h_u2.mat')['b1_2D']
    b20['g'][0]     = 1.*file['b2'][ics:ice,jcs:jce]                    # b_eff*sio.loadmat('input_b_h_u2.mat')['b2_2D']  
    # print('restart reading h10 and Oc from Topo_non_unigrid.mat')    
    # h10['g'][0]   = itopo*(sio.loadmat('Topo_non_unigrid.mat')['h_topo_smooth'][ics:ice,jcs:jce])/H0_dim
    # Oc['g'][0]    = sio.loadmat('Topo_non_unigrid.mat')['Ocean'][ics:ice,jcs:jce]
    # w1['g'][0]    = 0.*file['q1'][ics:ice,jcs:jce]#file['w1'][ics:ice,jcs:jce]    
    # w2['g'][0]    = 0.*file['q1'][ics:ice,jcs:jce]#file['w2'][ics:ice,jcs:jce]       
    oldt            =    file['t']
    t               = oldt[0] #0
#------------ end of restart initialization-------------#

save_data(comm, os.path.join(output_folder, 'debug_init_data.npz'),
          h1        = h1['g'][0], 
          h2        = h2['g'][0], 
          u1ph      = u1['g'][1], 
          u1th      = u1['g'][0], 
          u2ph      = u2['g'][1], 
          u2th      = u2['g'][0], 
          q1        = q1['g'][0], 
          b1        = b1['g'][0], 
          b2        = b2['g'][0], 
          q2        = q2['g'][0], 
          w1        = w1['g'][0], 
          w2        = w2['g'][0], 
          b10       = b10['g'][0], 
          b20       = b20['g'][0]
          # h10     = h10['g'][0]
          # Oc      = Oc['g'][0]
)

state_vector.pack(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)

# build matrices
P, M, L, LU = [], [], [], []

for m in range(m_start,m_end+1):
    Mm,Lm = eq.shallow_water(S,m,[g0, B3, B2, B1, H1, H2, Om, a, nu, kappa, Q01, Q02])
    M.append(Mm.astype(np.complex128))
    L.append(Lm.astype(np.complex128))
    P.append(0.*Mm.astype(np.complex128))
    # LU.append([None])

print0()
print0('Running --> Calculate RHS nonlinear terms from state_vector')

# calculate RHS nonlinear terms from state_vector
def nonlinear(state_vector,RHS):

    state_vector.unpack(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)

    Du1.layout      = 'c'
    Du2.layout      = 'c'
    Db1.layout      = 'c'
    Db2.layout      = 'c'    
    Dh1.layout      = 'c'
    Dh10.layout     = 'c'    
    Dh2.layout      = 'c'
          
    for m in range(m_start,m_end+1):
        md = m - m_start
        S.grad(m, 1,  u1['c'][md],  Du1['c'][md])
        S.grad(m, 1,  u2['c'][md],  Du2['c'][md])
        S.grad(m, 0,  b1['c'][md],  Db1['c'][md])
        S.grad(m, 0,  b2['c'][md],  Db2['c'][md])        
        S.grad(m, 0,  h1['c'][md],  Dh1['c'][md]) 
        S.grad(m, 0, h10['c'][md], Dh10['c'][md])         
        S.grad(m, 0,  h2['c'][md],  Dh2['c'][md])                      
                               
    u1_rhs.layout   = 'g'
    u2_rhs.layout   = 'g'
    h1_rhs.layout   = 'g'
    h2_rhs.layout   = 'g'
    uh1.layout      = 'g'
    uh2.layout      = 'g'   
    #uhf1.layout    = 'g'
    q1_rhs.layout   = 'g'
    q2_rhs.layout   = 'g'
    w1_rhs.layout   = 'g'
    w2_rhs.layout   = 'g'         
    uq1.layout      = 'g'
    uq2.layout      = 'g'
    uw1.layout      = 'g'
    uw2.layout      = 'g'  
    #uwf1.layout    = 'g'
    #uwf2.layout    = 'g'          
    #uqf1.layout    = 'g'
    #uqf2.layout    = 'g'     
    b1_rhs.layout   = 'g'
    b2_rhs.layout   = 'g'    
    # ub1.layout    = 'g'
    #ubf1.layout    = 'g'
    Ev.layout       = 'g'
    Oc.layout       = 'g'
    CC1.layout      = 'g'
    DD1.layout      = 'g'
    DDr.layout      = 'g'
    VV2.layout      = 'g'    
    VV1.layout      = 'g'   
    Prec1.layout    = 'g'  
    if Newtonian_cooling == 1:    
        Rad1.layout = 'g'
        Rad2.layout = 'g'
    RT1.layout      = 'g'
    RT2.layout      = 'g'
    #q_o_h1.layout  = 'g'       
    #q_o_h2.layout  = 'g'           
    #speed.layout   = 'g'
    R2.layout       = 'g'
    #Press1.layout  = 'g'
    #Press2.layout  = 'g'
    
    ## limit q1 to guarantee m_enthalpy>0 everywhere, q1 is humidity anomaly
    # q1['g'][0]          = np.minimum(q1['g'][0], (h1['g'][0]+H1)/ep1-1.0001*Q01) 
    # q1['g'][0]          = np.minimum(q1['g'][0], (h1['g'][0]+H1)/(ep1+gamma-0.1)-1.0001*Q01) 
    # if topography == 1:
        # q1['g'][0]      = np.minimum(q1['g'][0],  max_q) 
        # if i==0:
        #     q2['g'][0]= np.minimum(q2['g'][0],  max_q-Qs2) 
        # q2['g'][0]      = np.minimum(q2['g'][0],  max_q)  
        # if fast_run == 1:
            # q1['g'][0]  = np.minimum(q1['g'][0],  0.03) 
            # q2['g'][0]  = np.minimum(q2['g'][0],  0.03)          
        # if super_fast == 1:
            # q1['g'][0]  = np.minimum(q1['g'][0],  0.03) 
            # q2['g'][0]  = np.minimum(q2['g'][0],  0.03)                  
        # if uniform_DD == 0:    
            # DD1['g'][0] = np.minimum(DD1['g'][0], 0.002)
    # else:
        # q1['g'][0]      = np.minimum(q1['g'][0],  0.05)    
        # q2['g'][0]      = np.minimum(q2['g'][0],  0.05)  
        # if uniform_DD == 0:
            # DD1['g'][0] = np.minimum(DD1['g'][0], 0.002)
    b1_mean               = global_mean(comm, b1['g'][0])           
    T_s                   = T_s0*(B1+ b1_mean) 

    # -----------------------------------------------------------------------------
    # -------------- *** Calculation of variables for the coupler *** -------------
    # -----------------------------------------------------------------------------

    # Pressure in the lower and upper layer under weak gradient of potential temperature   
    # Press1['g'][0]    = p_level1 + g_real*p_level1*((h1['g'][0]+h2['g'][0]+h10['g'][0])*b1['g'][0])/a
    # Press2['g'][0]    = p_level2 + g_real*p_level2*( h1['g'][0]*b1['g'][0] +  (h10['g'][0]+h2['g'][0])*b2['g'][0])/a # h1b1+(h0+h2)b2
    # dPress2['g'][0]   = h1['g'][0]*Db1['g'][0]+Dh1['g'][0]*b1['g'][0]+0.5*h2['g'][0]*Db2['g'][0]+b2['g'][0]*(Dh10['g'][0]+h2['g'][0]) + h1['g'][0]*Db1['g'][1]+Dh1['g'][1]*b1['g'][0]+0.5*h2['g'][0]*Db2['g'][1]+b2['g'][0]*(Dh10['g'][1]+h2['g'][0])###-b2['g'][0]+h2['g'][0]
    f_mf                = q1['g'][0]/(global_amax(comm, q1['g'][0])+0.0000001) 
    R_real              = R_dry*(1.0 - f_mf) + R_moist*f_mf
    # gamma_trop        = 0.00649                #  lapse rate of the troposphere
    # T_bot             = T_s0*(b1['g'][0]+B1) - (g_real/R_real)*ma.log((1.0-delta1)/2.0) - gamma_trop * HORO#  HORO [m]
    # u_bot             = u_asr*u1['g'][1]*U_scale #[m/s]
    # v_bot             = u_asr*u1['g'][1]*U_scale #[m/s]
    # ----------- end of calculation of variables for the coupler -------------------


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #
    #                                                                    #
    # --------------- Calculation of Radiative Transfer ---------------- #
    #                                                                    #    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #        
    
    if Radiative_Transfer == 1:    
        stateT          = make_idealized_column(T_s)  # stateT['Tatm'] defines an idealized vertical structure (num_lev*1) of temperature 
        h2o             = climlab.radiation.water_vapor.ManabeWaterVapor(\
                          state             = stateT, 
                          relative_humidity = 0.8,
                          qStrat            = 5e-06)  #For dry case RH=0. and qStrat=0.                     
        # print('Temperature', stateT['Tatm'])
        # print('H2o', h2o.q)     # 1D columnar variation of specific humidity from O(10^-2) to O(10^-6)
        absorber_vmr = {'CO2'    :CO2ppmv/1e6,
                        'CH4'    :0.,
                        'N2O'    :0.,
                        'O2'     :0.,
                        'CFC11'  :0.,
                        'CFC12'  :0.,
                        'CFC22'  :0.,
                        'CCL4'   :0.,
                        'O3'     :0.}
        # RRTMG_LW radiation
        rad_LW          = climlab.radiation.RRTMG_LW(\
                          state               = stateT, 
                          specific_humidity   = h2o.q,        # h2o.q is 1D columnar variation of specific humidity from O(10^-2) to O(10^-6)
                          icld                = 2,            # Cloud overlap method, 0: Clear-sky only, 1: Random, 2: Maximum/random] 
                          return_spectral_olr = False,        # False: Just returns total OLR
                          ssac                = 0.,           # In-cloud single scattering albedo
                          asmc                = 0.,           # In-cloud asymmetry parameter, #0.6
                          fsfc                = 0.,           # In-cloud forward scattering fraction (delta function pointing forward "forward peaked scattering")
                          # ---- *** AEROSOLS *** ----- #
                          iaer                = 0,            # Aerosol option flag, 0: No aerosol, 6: ECMWF method: use six ECMWF aerosol types input aerosol optical depth at 0.55 microns for each aerosol type (ecaer), 10:Input aerosol optical properties: input total aerosol optical depth, single scattering albedo and asymmetry parameter (tauaer, ssaaer, asmaer) directly.
                          tauaer              = 0.,           # Aerosol optical depth (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
                          ssaaer              = 0.,           # Aerosol single scattering albedo (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
                          asmaer              = 0.,           # Aerosol asymmetry parameter (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
                          ecaer               = 0.,           # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions,  (ncol,nlay,naerec)] #  (non-delta scaled)   
                          absorber_vmr        = absorber_vmr)
        rad_LW.compute_diagnostics()
        rad_LW.OLR
        # print('OLR:      ', rad_LW.OLR)
        # print('OLR.shape:', rad_LW.OLR.shape)        
        # ------------------- *** Short Wave Radiative Transfer *** ---------------
        rad_SW          = climlab.radiation.RRTMG_SW(\
                          state             = stateT, 
                          specific_humidity   = h2o.q,        # h2o.q is 1D columnar variation of specific humidity from O(10^-2) to O(10^-6)
                          icld                = 2,            # Cloud overlap method, 0: Clear-sky only, 1: Random, 2: Maximum/random] 
                          ssac                = 0.,           # In-cloud single scattering albedo
                          asmc                = 0.,           # In-cloud asymmetry parameter, #0.6
                          fsfc                = 0.,           # In-cloud forward scattering fraction (delta function pointing forward "forward peaked scattering")
                          # ---- *** AEROSOLS *** ----- #
                          iaer                = 0,            # Aerosol option flag, 0: No aerosol, 6: ECMWF method: use six ECMWF aerosol types input aerosol optical depth at 0.55 microns for each aerosol type (ecaer), 10:Input aerosol optical properties: input total aerosol optical depth, single scattering albedo and asymmetry parameter (tauaer, ssaaer, asmaer) directly.
                          tauaer              = 0.,           # Aerosol optical depth (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
                          ssaaer              = 0.,           # Aerosol single scattering albedo (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
                          asmaer              = 0.,           # Aerosol asymmetry parameter (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
                          ecaer               = 0.,           # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions,  (ncol,nlay,naerec)] #  (non-delta scaled)   
                          absorber_vmr        = absorber_vmr)
        rad_SW.compute_diagnostics()
        #ncout.WriteVar('ASR', rad_SW.ASR, t*T_scale)      
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #
    #                                                                    #
    # ----------- End of Calculation of Radiative Transfer ------------- #
    #                                                                    #    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #    
        
    maxspeed            = global_amax(comm, np.abs(u1['g'][1])) + 0.000001 #global_amax(comm, speed['g'][0]) + 0.00001  
    minspeed            = global_amin(comm, np.abs(u1['g'][1])) 
    #maxspeed           = global_amax(comm, u1['g'][1]) + 0.000001 #global_amax(comm, speed['g'][0]) + 0.00001  
    #minspeed           = global_amin(comm, u1['g'][1])     
    b1_max              = global_amax(comm, b1['g'][0]) + 0.0000001
    #print('Maximum speed = ', maxspeed)
    #print('b1_min = '       , b1_min)
    #print('b1_max = '       , b1_max)
    # e_max             = global_amax(comm, Ev['g'][0])
    if Radiative_Transfer==1:      
        ti              = ma.floor(t*T_scale + t_init_insol)%365
        month_number    = ma.ceil((ti/365.0) * 12.0)
        #print('month_number',month_number)
        albedo_2Dt      = sio.loadmat('albedo_smooth.mat')['albedo_smooth'][ics:ice,jcs:jce,int(month_number)-1]
        albedo_crt      = 0.6
        
    # imposing external forcing
    '''                  
    tot_step   = 100
    step_size  = 5
    fb1_coeff  = 0.002#0.02
    b_coeff    = 0.8
    init       = 10 # The iteration number 1600 corresponds approximately to day 6.
    if i==1:
        b1['g'][0] = b_coeff*b1['g'][0] 
        b2['g'][0] = b_coeff*b2['g'][0] 
    if i in range(init, step_size*tot_step+init, step_size): # time span of external forcing for buoyancy anomaly
        if i==init:
            h1['g'][0] = h1['g'][0] - f1_ext 
            b1['g'][0] = b1['g'][0]+(1./tot_step)*fb1_coeff*f1_ext 
            print('Min of h1',global_amin(comm, h1['g'][0]))
            print('i = ',i)               
        else:
            b1['g'][0] += (1./tot_step)*fb1_coeff*f1_ext
            print('Max of b1',global_amax(comm, b1['g'][0]))
            print('Min of b1',global_amin(comm, b1['g'][0]))
            print('i = ',i)
    '''
        
    if itopo == 1.:
        # TODO: Surface evaporation over the lands needs to be modified.

        alpha1v           = 1.2       # 1.2, higher: higher peak of CLWC in off-equatorial region with high amplitude of velocity OR more separation of CLWC and q between low and high velocity.
        alpha2v           = 0.7       # 0.7, higher: extends latitudinal propagation of water vapor and CLWC & more dependency on high velocity
        u_norm            = (np.abs(u1['g'][1])-minspeed)/(maxspeed-minspeed)        
        #e_v              = np.heaviside(((b1['g'][0]+B1)*T_s0)-T0,0.) * np.exp((u_norm**alpha1v)/alpha2v) #+ np.exp((np.abs(u1['g'][0]) /maxspeed)**alpha1v/alpha2v)   
        e_v               = np.exp((u_norm**alpha1v)/alpha2v)         
        if t*T_scale<spin_up_days:          
            coeff_b       = 0.4     # 
            coeff_v       = 0.04     
            coeff_sf_ev   = 0.5     # 0.4  scaling factor correction of the Clausius-Clapeyron equation
            coeff_sf_qs   = 0.5     # 0.4  scaling factor correction of the Clausius-Clapeyron equation              
        else:
            coeff_b       = 0.055   # 0.05
            coeff_v       = 0.0005  # 0.0005
            coeff_sf_ev   = 0.4     # 0.4  scaling factor correction of the Clausius-Clapeyron equation
            coeff_sf_qs   = 0.4     # 0.4  scaling factor correction of the Clausius-Clapeyron equation        

        free_conv         = 0.0     # 0.00001 free convection over the ocean, i.e. evaporation when there is no wind, no heating, hogher leads to sharper at Eq.
        exp_arg           = (-L_v / R_moist) * (np.clip(1.0 / np.abs((b1['g'][0] + B1 + 1e-6) * T_s0), a_min=None, a_max=1e10) ** coeff_sf_ev - 1.0 / (T0 ** coeff_sf_ev))
        
        exp_arg           = np.where(exp_arg < -13.81, np.NINF, exp_arg)  # np.log(1e-6) = -13.81  or  np.exp(-13.81) = 1e-6  AND  np.exp(np.NINF) = 0.0
        exp_arg           = np.where(exp_arg >  59.86,   59.86, exp_arg)  # np.log(1e26) =  59.86  or  np.exp( 59.86) = 1e26  
        
        e_b               = np.exp(exp_arg)     # before --> np.where(np.abs(exp_arg) < 1e-15, 0.0, np.exp(exp_arg))
                
        e_b              *= np.heaviside(((b1['g'][0] + B1) * T_s0) - T0, 0.) * np.heaviside(b1['g'][0] + B1, 0.)
        #e_b              = np.heaviside(((b1['g'][0] + B1) * T_s0) - T0, 0.) * np.exp((-L_v / R_moist) * (np.clip(1.0 / np.abs((b1['g'][0] + B1 + 1e-6) * T_s0), a_min=None, a_max=1e10) ** coeff_sf_ev - 1.0 / (T0 ** coeff_sf_ev))) * np.heaviside(b1['g'][0] + B1, 0.)
        #e_b              = np.heaviside(((b1['g'][0]+B1)*T_s0)-T0,0.) * np.exp((-L_v/R_moist) * (1.0/(np.abs((b1['g'][0]+B1+0.00000001)*T_s0)**coeff_sf_ev) - 1.0/(T0**coeff_sf_ev))) * np.heaviside(b1['g'][0]+B1,0.)           
        e_v_max           = global_amax(comm, e_v) +0.00000001
        e_v_min           = global_amin(comm, e_v) 
        e_b_max           = global_amax(comm, e_b) +0.00000001
        e_b_min           = global_amin(comm, e_b)         
        e_vv              = np.heaviside(albedo_crt-albedo_2Dt, 0.0) *((e_v-e_v_min)/(e_v_max-e_v_min))         
        e_bb              = (e_b-e_b_min)/(e_b_max-e_b_min) 
        Qs_limit          = 1.2* Qs1 #1.3#
        #print('Qs_limit = ', global_amax(comm, Qs_limit))
        
        '''
        aa1                 = q1['g'][0]
        aa2                 = q2['g'][0]        
        aa1[Qs_limit < aa1] = Qs_limit[Qs_limit < aa1] 
        aa2[Qs_limit < aa2] = Qs_limit[Qs_limit < aa2]  
        q1['g'][0]          = aa1
        q2['g'][0]          = aa2
        '''
        
        q1['g'][0]    = np.minimum(q1['g'][0],  Qs_limit)
        q2['g'][0]    = np.minimum(q2['g'][0],  Qs_limit)
        
        #print('size of Qs1 = '       , Qs1.shape)
        #print('e_vv_max = '       , global_amax(comm, e_vv))
        # e_bv        = (e_bb*e_vv) 

        # Ev['g'][0]  = (i_mc*Oc['g'][0]*np.heaviside(f_r*max_q-q1['g'][0]-Q01,0.) *np.minimum(f_r*max_q-q1['g'][0]-Q01, f_r*max_q)*e_bv)+ (i_mc*np.heaviside(f_r*max_q-q1['g'][0]-Q01,0.)* (1.-Oc['g'][0])  * np.minimum(f_r*max_q-q1['g'][0]-Q01, f_r*max_q)*e_bv/10.0)  
        
        #Ev['g'][0]   = np.heaviside(((Qs_limit+Qs1)/2.0)-q1['g'][0]-Q01,0.)*(coeff_b*((i_mc*Oc['g'][0]*(e_bb+free_conv))+ (i_mc* (1.-Oc['g'][0])*e_bb/(15.0))) + coeff_v*((i_mc*Oc['g'][0]*e_vv)+ (i_mc* (1.-Oc['g'][0])*e_vv/2.2)))  #1.0   
        Ev['g'][0]    = np.heaviside(Qs1-q1['g'][0]-Q01,0.)  *  (coeff_b*((i_mc*Oc['g'][0]*(e_bb+ (np.heaviside(albedo_crt-albedo_2Dt, 0.0)*free_conv)))) + coeff_v*((i_mc*Oc['g'][0]*e_vv)+ (i_mc* (1.-Oc['g'][0])*e_vv/2.5)))  #1.7                                                                          
        CC1['g'][0]   = i_mc*np.heaviside(q1['g'][0]+Q01-Qs1,0.)*(q1['g'][0]+Q01-Qs1)*t_p_inv   
        #CC1['g'][0]  = i_mc*np.heaviside(q1['g'][0]+Q01-Qs1,0.)*np.minimum((q1['g'][0]+Q01-Qs1), 2.0*Qs1)*t_p_inv               
    else:
        q1['g'][0]    = np.minimum(q1['g'][0],  0.05)    
        q2['g'][0]    = np.minimum(q2['g'][0],  0.05)  
        if uniform_DD == 0:
            DD1['g'][0] = np.minimum(DD1['g'][0], 0.002)    
        # CC1['g'][0] = i_mc*np.heaviside(q1['g'][0]+Q01-H1*Qs1,0.)*np.minimum((q1['g'][0]+Q01-H1*Qs1), 1.5*H1*Qs1)*t_p_inv #0.015, 0.03, normal case   
        CC1['g'][0]   = i_mc*np.heaviside(q1['g'][0]+Q01-Qs1,0.)*np.minimum((q1['g'][0]+Q01-Qs1), 0.05)*t_p_inv 

        Ev['g'][0]    = i_mc*np.heaviside(Qs1-q1['g'][0]-Q01,0.)*(Qs1-q1['g'][0]-Q01)*np.minimum(Qs1-q1['g'][0]-Q01, 0.015)*t_e_inv*(np.abs(u1['g'][1])/maxspeed)
        # CC1['g'][0] = i_mc*np.heaviside(q1['g'][0]+Q01-Qs1,0.)*np.minimum((q1['g'][0]+Q01-Qs1), 0.05)*t_p_inv 
        # Ev can be proportional to: (b-b_e)/b0
    
    
    heavQ = np.heaviside(Qs1-(q1['g'][0]+Q01), 0.)
    if uniform_DD == 1:
        nonzeroDDnum   = global_count_nonzero(comm, h1)
        DD1['g'][0]    = heavQ*global_sum(comm, CC1['g'][0]-Ev['g'][0])/nonzeroDDnum # it can be parametrized with respect to b2, once initial b2 is not zero
        DDr['g'][0]    = heavQ*global_sum(comm, CC1['g'][0])/nonzeroDDnum
        q2['g'][0]     = q2['g'][0] + ((global_sum(comm, Ev['g'][0]))/(Nxx*Nyy)) # This balances the loss of water vapor in the upper layer
        # *** Bias correction of q1 ***
        EVCCDDsum      = global_sum(comm, (Ev['g'][0]-CC1['g'][0]+DD1['g'][0]))
        q1sum          = global_sum(comm, q1['g'][0])
        bias_q1        = (EVCCDDsum-q1sum)/(Nxx*Nyy)
        q1['g'][0]     = q1['g'][0] + iDD*bias_q1    
    else: # nonuniform downdraft    ***
        #DD1['g'][0]    = heavQ* np.minimum((Qs1-(q1['g'][0]+Q01)), 0.995*Qs1)*t_p_inv#0.002
        # Check if each element of the 2D variable q1['g'][0] falls within the specified range
        lower_limit    = 0.65*Qs1
        upper_limit    = 0.999*Qs1
        DD1['g'][0]    = heavQ*np.logical_and(q1['g'][0] >= lower_limit, q1['g'][0] <= upper_limit)
        # np.logical_and() is a NumPy function that computes the element-wise logical AND operation between two arrays or between an array and a scalar.
        # True =1, False = 0
        DD1sum         = global_sum(comm, DD1['g'][0])+0.0000001
        #print('Coefficient of D1, old coeff', Coeff_D1, Coeff_D1_bis)     
        Coeff_D1       = global_sum(comm, (CC1['g'][0]-Ev['g'][0]))/(DD1sum+0.0000001) # this coefficient closes the cycle of water vapor, so helps to conservation of q1       
        Coeff_Dr       = global_sum(comm, CC1['g'][0])/DD1sum
        DDr['g'][0]    = Coeff_Dr*DD1['g'][0]   # without closed water cycle
        DD1['g'][0]    = iDD*DDr['g'][0]  #iDD*Coeff_D1*DD1['g'][0]
        q2['g'][0]     = q2['g'][0] + iDD*global_sum(comm, Ev['g'][0])/(Nxx*Nyy) # This balances the loss of water vapor in the upper layer, which is just a passive tracer       

        # *** Bias correction of q1 ***
        # EVCCDDsum    = global_sum(comm, (Ev['g'][0]-CC1['g'][0]+iDD*DD1['g'][0]))
        EVCCDDsum      = global_sum(comm, (Ev['g'][0]-CC1['g'][0]))
        q1sum          = global_sum(comm, q1['g'][0])

        bias_q1        = (EVCCDDsum-q1sum)/(1.0*Nxx*Nyy)
        # bias_E1      = EVCCDDsum/(0.6*Nxx*Nyy)
        # Ev['g'][0]   = Ev['g'][0] - bias_E1*Oc['g'][0]
        # q1['g'][0]   = q1['g'][0] + bias_q1#*Oc['g'][0]
       
       
    '''              
    if uniform_DD == 1:
        nonzeroDDnum   = np.count_nonzero(np.heaviside(Qs1-(q1['g'][0]+Q01),0.))
        DD1['g'][0]    = np.heaviside(Qs1-(q1['g'][0]+Q01),0.)*(CC1['g'][0]-Ev['g'][0]).sum()/nonzeroDDnum # it can be parametrized with respect to b2, once initial b2 is not zero
        DDr['g'][0]    = np.heaviside(Qs1-(q1['g'][0]+Q01),0.)*(CC1['g'][0]).sum()/nonzeroDDnum  
        q2['g'][0]     = q2['g'][0] + (((Ev['g'][0]).sum())/(lamda.size*phi.size)) # This balances the loss of water vapor in the upper layer
        # *** Bias correction of q1 ***
        if q1['g'][0].sum()<((Ev['g'][0]-CC1['g'][0]+DD1['g'][0]).sum()):
            bias_q1    = (((Ev['g'][0]-CC1['g'][0]+DD1['g'][0]).sum())-q1['g'][0].sum())/(lamda.size*phi.size)
            q1['g'][0] = q1['g'][0] + bias_q1
        if q1['g'][0].sum()>((Ev['g'][0]-CC1['g'][0]+DD1['g'][0]).sum()):
            bias_q1    = (-((Ev['g'][0]-CC1['g'][0]+DD1['g'][0]).sum())+q1['g'][0].sum())/(lamda.size*phi.size)         
            q1['g'][0] = q1['g'][0] - bias_q1       
    else: # nonuniform downdraft   ***
        DD1['g'][0]    = np.heaviside(Qs1-(q1['g'][0]+Q01),0.)* np.minimum((Qs1-(q1['g'][0]+Q01)), 0.002)*t_p_inv  
        #print('Coefficient of D1, old coeff', Coeff_D1, Coeff_D1_bis)     
        #Coeff_D1      = (CC1['g'][0]-Ev['g'][0]).sum()/(DD1['g'][0]).sum() # this coefficient closes the cycle of water vapor, so helps to conservation of q1
        Coeff_D1       = global_sum(comm, (CC1['g'][0]-Ev['g'][0]))/(global_sum(comm, DD1['g'][0])+0.0000001)# this coefficient closes the cycle of water vapor, so helps to conservation of q1       
        Coeff_Dr       = (CC1['g'][0]).sum()/(DD1['g'][0]).sum() 
        DDr['g'][0]    = Coeff_Dr*DD1['g'][0]   # without closed water cycle
        DD1['g'][0]    = Coeff_D1*DD1['g'][0]
        q2['g'][0]     = q2['g'][0] + global_sum(comm, Ev['g'][0])/(Nxx*Nyy) # This balances the loss of water vapor in the upper layer       
        #q2['g'][0]    = q2['g'][0] + (((Ev['g'][0]).sum())/(lamda.size*phi.size)) # This balances the loss of water vapor in the upper layer
        # *** Bias correction of q1 ***
        if global_sum(comm, q1['g'][0])<global_sum(comm, (Ev['g'][0]-CC1['g'][0]+DD1['g'][0])):
            bias_q1    = (global_sum(comm, Ev['g'][0]-CC1['g'][0]+DD1['g'][0]) -  global_sum(comm,q1['g'][0]))  /(Nxx*Nyy)
            q1['g'][0] = q1['g'][0] + bias_q1
        if global_sum(comm, q1['g'][0])>global_sum(comm, (Ev['g'][0]-CC1['g'][0]+DD1['g'][0])):
            bias_q1    = (-global_sum(comm, Ev['g'][0]-CC1['g'][0]+DD1['g'][0])+global_sum(comm,q1['g'][0]))/(Nxx*Nyy)        
            q1['g'][0] = q1['g'][0] - bias_q1        
    '''                            
    
    VV2['g'][0]     = ivap*w2['g'][0]*t_v_inv    
    VV1['g'][0]     = ivap*w1['g'][0]*t_v_inv  
    
    Prec1['g'][0]   = np.heaviside(w1['g'][0]-wcritical,0.)*(w1['g'][0]-wcritical)*0.8*t_p_inv                                                               
    # print('Total q: ',(q1['g'][0]).sum(),'  Es:', (Ev['g'][0]).sum(),'  C1:', (CC1['g'][0]).sum(),'  D1:', (DD1['g'][0]).sum())
    
    R2['g'][0]        = ( ep1*(CC1['g'][0]-iDD*DD1['g'][0])*(u1['g'][0]-u2['g'][0]) + 1e-10 ) / ( (H2+h2['g'][0]*(b2['g'][0]+B2)) + 1e-10 )
        
    # Stokes drag that can be implemented in the upper or lower layer due to condensation in the lower one
    R2['g'][1]        = ( ep1*(CC1['g'][0]-iDD*DD1['g'][0])*(u1['g'][1]-u2['g'][1]) + 1e-10 ) / ( (H2+h2['g'][0]*(b2['g'][0]+B2)) + 1e-10 )
    
    # [(u1-u2)/(h2)]*W(velocity at the interface due to heat transfer)
    if Newtonian_cooling == 1:
       Rad1['g'][0] = -1.0*iradN*(((h1['g'][0]+H1)*(b1['g'][0]+B1))-H1*(B1+b10['g'][0]))*t_r_inv-iradh*h1['g'][0]*t_r_inv 
       Rad2['g'][0] = -1.0*iradN*(((h2['g'][0]+H2)*(b2['g'][0]+B2))-H2*(B2+b20['g'][0]))*t_r_inv-iradh*h2['g'][0]*t_r_inv  
    if   summer_solstice == 1:
        insol2D     = iRT*np.tile(insol[:,(177)], (lons.size, 1)) 
    elif winter_solstice == 1:
        insol2D     = iRT*np.tile(insol[:,(353)], (lons.size, 1))                   
    else:
        insol2D     = iRT*np.tile(insol[:,ma.floor(t*T_scale + t_init_insol)%365], (lons.size, 1)) 
        insol2D_ref = iRT*np.tile(insol[:,ma.floor(0.0       + t_init_insol)%365], (lons.size, 1))
        # Note that when selecting 'insol2D_ref,' the external forcing of radiative transfer will be based on the insolation anomaly relative to the chosen day
    # print('insol2D shape', insol2D.shape)
    # ncout.WriteVar('insolation',  insol2D[0,:], t*T_scale)     # doesn't work
    Q_max_global    = 565. # maximum global value of the daily insol (according to the Climlab data) which is not essentially equal to max(insol) during parallel running
    
    # Q_max_global  = global_amax(comm, insol2D)
    # print('*** Maximum Global INsolation ***', Q_max_global)
    
    
    if Radiative_Transfer==1:
       RT2Di                 = iRT*insol2D*(1.0+(rad_SW.ASR-rad_LW.OLR)/Q_max_global)      # the OLR can never reach the input insolation
       if Insolation_Exclusive:
          #ampl_RT2D         = global_amax(comm, RT2Di) - global_amin(comm, RT2Di)
          coeff_RT1          = 1.0# 100.0# 1.275#1.7#2.0#1.5#  normal: 1.5, Eq_warming: 2, <2.0
          coeff_RT2          = 1.0# 100.0# 1.0# 1.275#1.7#2.0#1.5#
          #anomaly_adjustment = 0.9*ampl_RT2D #0.8<...<1.0   
          #RT2D              = ampl_RT2D*(1.-albedo_2Dt)*(RT2Di-anomaly_adjustment)  
          RT2D               = (1.0-albedo_2Dt)*RT2Di            
       elif Anomaly_Mean_Critical==1:          # 1: Indicates a radiative transfer anomaly significantly deviating from the mean value.  
          ampl_RT2D          = global_amax(comm, RT2Di) - global_amin(comm, RT2Di)
          coeff_RT1          = 0.07#0.10
          coeff_RT2          = 0.07#0.10  
          thermalDissipation = 2.3*ampl_RT2D      
          RT2D_ref1          = 0.9*global_mean(comm, insol2D_ref) #0.1*global_mean(comm, RT2Di) #
          RT2D               = ampl_RT2D*(1.0-albedo_2Dt)*(np.heaviside(RT2Di-RT2D_ref1,0.)*(RT2Di-RT2D_ref1+thermalDissipation) + np.heaviside(-RT2Di+RT2D_ref1,0.)*(RT2Di-RT2D_ref1-thermalDissipation))
          #print('global ampl. of RT2D',ampl_RT2D)    
          #print('global mean of RT2D',global_mean(comm, RT2Di))   
       '''   
       elif Anomaly_Ref_Critical==1:
          coeff_RT1          = 0.05# <1.
          coeff_RT2          = 0.05   
          thermalDissipation = 1.9*ampl_RT2D #250.0*nu 
          RT2D_ref2          = iRT*insol2D_ref*(1.0+(rad_SW.ASR-rad_LW.OLR)/Q_max_global)
          RT2D               = ampl_RT2D*((1.0-albedo_2Dt)*np.heaviside(RT2Di-RT2D_ref2,0.)*(RT2Di-RT2D_ref2+thermalDissipation) + np.heaviside(-RT2Di+RT2D_ref2,0.)*(RT2Di-RT2D_ref2-thermalDissipation))  
      '''            
          #RT2D              = ampl_RT2D*(RT2Di-RT2D_ref2) #  Radiative transfer anomaly with respect to the reference level (day)
          #print('global max. of RT2D', global_amax(comm, RT2D)) 
          #print('albedo_2Dt shape', albedo_2Dt.shape)  
          #print('RT2Di shape', RT2Di.shape)  
          #print('month_number = ', month_number)
    else:
        RT2D                 = iRT*insol2D           
    # print('RT2D', RT2D) 
      
    rho_atm1        = 1.0
    rho_atm2        = 0.8                       # Average density of air  (kg/m^3) with respect to the vertical altitude above sea level on H2 (requires modification).
    # H_atm         = 10000.0
    C_dry           = 1005.0                    # Specific  heat capacity (J/kg/K) of dry air at constant pressure.
    C_vap           = 1875.0                    # Specific  heat capacity (J/kg/K) of water vapor at constant pressure. 
    C_eff           = f_mf*C_vap+(1.0- f_mf)*C_dry
    RT1['g'][0]     = coeff_RT1*B1*(RT2D/(C_eff*delta1*rho_atm1))*(1.0/T_s0)  # delta_b = delta T/thta_s: buoyancy increase 
    RT2['g'][0]     = coeff_RT2*B2*(RT2D/(C_eff*delta2*rho_atm2))*(1.0/T_s0)  #    
    #if i in range(1, 25*10+1, 25): 
        #print('Max of RT1' ,global_amax(comm, RT1['g'][0]))
        #print('i = ',i)                                  
    #print('delta T Land1 {:.10f}'.format(np.max(deltaT_Lnd1)))
    # deltaT_O_Ts = (1./Co_Cl)*(1./(T_s))*RT2D  # daily T increase over the ocean                                     
    # RT1['g'][0] = coeff_RT*( (h1['g'][0]+H1)*(B1+(Oc['g'][0]*deltaT_Ocn1+(1.-Oc['g'][0])*deltaT_Lnd1)) - H1*B1)#         
    # RT2['g'][0] = coeff_RT*( (h2['g'][0]+H2)*(B2+deltaT_Lnd2) - H2*B2) # density of water vapor in the upper layer is set similar to land         
                 
    # RT1['g'][0] = coeff_RT*(1.0/g_real)*( (h1['g'][0]+H1)*(B1+(Oc['g'][0]*deltaT_O_Ts+Co_Cl*(1.-Oc['g'][0])*deltaT_O_Ts)) - H1*B1)#  
    # RT2['g'][0] = coeff_RT*(1.0/g_real)*( (h2['g'][0]+H2)*(B2+(Oc['g'][0]*deltaT_O_Ts+Co_Cl*(1.-Oc['g'][0])*deltaT_O_Ts)) - H2*B2)                         
    # print('max RT1',np.max(RT1['g'][0]))   # order= 10                  

    # deltaT_O_Ts = (1./T_s)*RT2D  # nondimensional form
        # print('max deltaT_O_Ts',np.max(deltaT_O_Ts))     
              
    # RT1['g'][0] = coeff_RT*( (h1['g'][0]+H1)*(B1+deltaT_Ocn1) - H1*B1)
    # RT2['g'][0] = coeff_RT*( (h2['g'][0]+H2)*(B2+deltaT_Lnd2) - H2*B2)        

    # RT1['g'][0] = coeff_RT*(1.0/g_real)*(   (h1['g'][0]+H1)*(B1+deltaT_O_Ts)  - H1*B1   )# for daily insolation      
    # RT2['g'][0] = coeff_RT*(1.0/g_real)*(   (h2['g'][0]+H2)*(B2+deltaT_O_Ts)  - H2*B2   )# for daily insolation      
                
    # print('RT1 shape:', RT1['g'][0].shape, 't floor', ma.floor(t/T_scale))
    # print('max RT1', np.max(RT1['g'][0]))
    # print('min RT',  np.min(RT1['g'][0]))
    #q_o_h1['g'][0]  = iq_o_h*((q1['g'][0]+Q01)/(h1['g'][0]+H1))
    #q_o_h2['g'][0]  = iq_o_h*((q2['g'][0]+Q02)/(h2['g'][0]+H2))
    
    
    
    # Adding delimiters for numerical stability and hyperbolicity is recommended.   
    h1['g'][0] = np.minimum(0.85 * delta1, h1['g'][0])  # Apply lower limit
    h1['g'][0] = np.maximum(-0.85* delta1, h1['g'][0])  # Apply upper limit 

    h2['g'][0] = np.minimum(0.85 * delta2, h2['g'][0])  # Apply lower limit
    h2['g'][0] = np.maximum(-0.85* delta2, h2['g'][0])  # Apply upper limit         
        
# ------------- *** Start of adding forces on the RHS of equations for three configurations *** -----------------#
    e_s_f      = global_sum(comm, CC1['g'][0])/(global_sum(comm, Ev['g'][0])+0.0000001)
    # evaporation_scaling_factor: A scaling factor applied to the surface/source evaporation rate to adjust its contribution to the overall moisture budget in the model.
    if smooth_run == 1:
        u1_rhs['g'][0]   = - (u1['g'][0]*Du1['g'][0] + u1['g'][1]*Du1['g'][2])/a - t_drag_inv*u1['g'][0] - 0.5*(h1['g'][0]*Db1['g'][0])/a-(b1['g'][0]*(Dh1['g'][0]+itopo*Dh10['g'][0]))/a - (b1['g'][0]*Dh2['g'][0])/a -  step*(u1['g'][0]-u10['g'][0])*t_r_u1
        u1_rhs['g'][1]   = - (u1['g'][0]*Du1['g'][1] + u1['g'][1]*Du1['g'][3])/a - t_drag_inv*u1['g'][1] - 0.5*(h1['g'][0]*Db1['g'][1])/a-(b1['g'][0]*Dh1['g'][1]+itopo*Dh10['g'][1])/a - (b1['g'][0]*Dh2['g'][1])/a  - step*(u1['g'][1]-u10['g'][1])*t_r_u1
        u2_rhs['g'][0]   = R2['g'][0] - (u2['g'][0]*Du2['g'][0] + u2['g'][1]*Du2['g'][2])/a - t_drag_inv*u2['g'][0] - 0.5*(h2['g'][0]*Db2['g'][0])/a-(b2['g'][0]*(Dh2['g'][0]+itopo*Dh10['g'][0]))/a - (h1['g'][0]*Db1['g'][0])/a - (b1['g'][0]*Dh1['g'][0])/a -  step*(u2['g'][0]-u20['g'][0])*t_r_u2
        u2_rhs['g'][1]   = R2['g'][1]  - (u2['g'][0]*Du2['g'][1] + u2['g'][1]*Du2['g'][3])/a - t_drag_inv*u2['g'][1] - 0.5*(h2['g'][0]*Db2['g'][1])/a-(b2['g'][0]*(Dh2['g'][1]+itopo*Dh10['g'][1]))/a - (h1['g'][0]*Db1['g'][1])/a - (b1['g'][0]*Dh1['g'][1])/a - step*(u2['g'][1]-u20['g'][1])*t_r_u2
        if Newtonian_cooling == 1:
           h1_rhs['g'][0]   = (1.0/(b1['g'][0]+B1))  * (ep1*(-CC1['g'][0]+iDD*DD1['g'][0]) + ((gamma_NC1-1.0)*Rad1['g'][0]) + ((gamma_RT1-1.0)*RT1['g'][0]))   
           h2_rhs['g'][0]   = (1.0/(b2['g'][0]+B2))  * (ep1*( CC1['g'][0]-iDD*DD1['g'][0]) + ((gamma_NC2-1.0)*Rad2['g'][0]) + ((gamma_RT2-1.0)*RT2['g'][0]))  
           b1_rhs['g'][0]   = - (u1['g'][0]*Db1['g'][0] + u1['g'][1]*Db1['g'][1])/a - step*(b1['g'][0]-b10['g'][0])*t_r_b1   +  step*(( CC1['g'][0]-iDD*e_s_f*Ev['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt   +   step*RT1['g'][0]/(h1['g'][0]+H1)/tau_bu*dt  + step*((Rad1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt                       
           b2_rhs['g'][0]   = - (u2['g'][0]*Db2['g'][0] + u2['g'][1]*Db2['g'][1])/a - step*(b2['g'][0]-b20['g'][0])*t_r_b2   +  step*((-CC1['g'][0]+iDD*DD1['g'][0]))/(h2['g'][0]+H2)/tau_bu*dt  +   step*RT2['g'][0]/(h2['g'][0]+H2)/tau_bu*dt + step*(Rad2['g'][0])/(h2['g'][0]+H2)/tau_bu*dt        
        else:       
           h1_rhs['g'][0]   = (1.0/(b1['g'][0]+B1))  * (ep1*(-CC1['g'][0]+iDD*DD1['g'][0])+((gamma_RT1-1.0)*RT1['g'][0]))          
           h2_rhs['g'][0]   = (1.0/(b2['g'][0]+B2))  * (ep1*( CC1['g'][0]-iDD*DD1['g'][0])+((gamma_RT2-1.0)*RT2['g'][0]))
           b1_rhs['g'][0]   = - (u1['g'][0]*Db1['g'][0] + u1['g'][1]*Db1['g'][1])/a - step*(b1['g'][0]-b10['g'][0])*t_r_b1   +  step*(( CC1['g'][0]-iDD*e_s_f*Ev['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt   +   step*RT1['g'][0]/(h1['g'][0]+H1)/tau_bu*dt                          
           b2_rhs['g'][0]   = - (u2['g'][0]*Db2['g'][0] + u2['g'][1]*Db2['g'][1])/a - step*(b2['g'][0]-b20['g'][0])*t_r_b2   +  step*((-CC1['g'][0]+iDD*DD1['g'][0]))/(h2['g'][0]+H2)  +   step*RT2['g'][0]/(h2['g'][0]+H2) 
           
           
           
        w1_rhs['g'][0]   = ivap*((1.0-upward)*CC1['g'][0] - Prec1['g'][0])#2nd_term=precipitation
        w2_rhs['g'][0]   = ivap*(upward*CC1['g'][0] - VV2['g'][0])        
        q2_rhs['g'][0]   = CC1['g'][0] - iDD*DD1['g'][0] # ivap*VV2['g'][0]#           #\nabla.{(uq2)} is added seperately 
        q1_rhs['g'][0]   = Ev['g'][0]  - CC1['g'][0] #+ iDD*DD1['g'][0]    #\nabla.{(uq1)} is added seperately 
                
        #b1_rhs['g'][0]   = - (u1['g'][0]*Db1['g'][0] + u1['g'][1]*Db1['g'][1])/a - step*(b1['g'][0]-b10['g'][0])*t_r_b1   +  step*(( CC1['g'][0]-iDD*DD1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt   +   step*RT1['g'][0]/(h1['g'][0]+H1)/tau_bu*dt #+ step*((Rad1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt                      
        #b2_rhs['g'][0]   = - (u2['g'][0]*Db2['g'][0] + u2['g'][1]*Db2['g'][1])/a - step*(b2['g'][0]-b20['g'][0])*t_r_b2   +  step*((-CC1['g'][0]+iDD*DD1['g'][0]))/(h2['g'][0]+H2)/tau_bu*dt  +   step*RT2['g'][0]/(h2['g'][0]+H2)/tau_bu*dt # +   step*(Rad2['g'][0])/(h2['g'][0]+H2)/tau_bu*dt        
         
        #b1_rhs['g'][0]   = - (u1['g'][0]*Db1['g'][0] + u1['g'][1]*Db1['g'][1])/a - step*(b1['g'][0]-b10['g'][0])*t_r_b1 + step*((1.0-(q_o_h1['g'][0]*ep1)/(B1+b1['g'][0]+q_o_h1['g'][0]))*(CC1['g'][0]-iDD*DD1['g'][0]) )/(h1['g'][0]+H1)/tau_bu*dt  +step*(    (  ((q_o_h1['g'][0]/(B1+b1['g'][0]+q_o_h1['g'][0])))  * ((gamma_NC1-1.0)*Rad1['g'][0]+(gamma_RT1-1.0)*RT1['g'][0])  )  +    (Rad1['g'][0]+RT1['g'][0])    )/(h1['g'][0]+H1)/tau_bu*dt                      
        #b2_rhs['g'][0]   = - (u2['g'][0]*Db2['g'][0] + u2['g'][1]*Db2['g'][1])/a - step*(b2['g'][0]-b20['g'][0])*t_r_b2 + step*((1.0-q_o_h2['g'][0]*ep1)*(-CC1['g'][0] + iDD*DD1['g'][0]))/(h2['g'][0]+H2)/tau_bu*dt     +step*(    (  ((q_o_h2['g'][0]/(B2+b2['g'][0]+q_o_h2['g'][0])))  * ((gamma_NC2-1.0)*Rad2['g'][0]+(gamma_RT2-1.0)*RT2['g'][0])  )  +    (Rad2['g'][0]+RT2['g'][0])    )/(h2['g'][0]+H2)/tau_bu*dt
        # print('t_r_b1', t_r_b1) 
        # print('max b10', np.max(b10['g'][0]))

              
    if fast_run == 1:
        u1_rhs['g'][0]   = - (u1['g'][0]*Du1['g'][0] + u1['g'][1]*Du1['g'][2])/a - t_drag_inv*u1['g'][0] - 0.5*(h1['g'][0]*Db1['g'][0])/a-(b1['g'][0]*(Dh1['g'][0]+itopo*Dh10['g'][0]))/a - (b1['g'][0]*Dh2['g'][0])/a -  step*(u1['g'][0]-u10['g'][0])*t_r_u1
        u1_rhs['g'][1]   = - (u1['g'][0]*Du1['g'][1] + u1['g'][1]*Du1['g'][3])/a - t_drag_inv*u1['g'][1] - 0.5*(h1['g'][0]*Db1['g'][1])/a-(b1['g'][0]*Dh1['g'][1]+itopo*Dh10['g'][1])/a - (b1['g'][0]*Dh2['g'][1])/a  - step*(u1['g'][1]-u10['g'][1])*t_r_u1
        u2_rhs['g'][0]   = R2['g'][0] - (u2['g'][0]*Du2['g'][0] + u2['g'][1]*Du2['g'][2])/a - t_drag_inv*u2['g'][0] - 0.5*(h2['g'][0]*Db2['g'][0])/a-(b2['g'][0]*(Dh2['g'][0]+itopo*Dh10['g'][0]))/a - (h1['g'][0]*Db1['g'][0])/a - (b1['g'][0]*Dh1['g'][0])/a -  step*(u2['g'][0]-u20['g'][0])*t_r_u2
        u2_rhs['g'][1]   = R2['g'][1]  - (u2['g'][0]*Du2['g'][1] + u2['g'][1]*Du2['g'][3])/a - t_drag_inv*u2['g'][1] - 0.5*(h2['g'][0]*Db2['g'][1])/a-(b2['g'][0]*(Dh2['g'][1]+itopo*Dh10['g'][1]))/a - (h1['g'][0]*Db1['g'][1])/a - (b1['g'][0]*Dh1['g'][1])/a - step*(u2['g'][1]-u20['g'][1])*t_r_u2 
        h1_rhs['g'][0]   = (1.0/(b1['g'][0]+B1))  * (ep1*(-CC1['g'][0]+iDD*DD1['g'][0])+((gamma_RT1-1.0)*RT1['g'][0])) #+((gamma_NC1-1.0)*Rad1['g'][0])           
        h2_rhs['g'][0]   = (1.0/(b2['g'][0]+B2))  * (ep1*( CC1['g'][0]-iDD*DD1['g'][0])+((gamma_RT2-1.0)*RT2['g'][0]))#+((gamma_NC2-1.0)*Rad2['g'][0])
        w1_rhs['g'][0]   = ivap*((1.0-upward)*CC1['g'][0] - Prec1['g'][0])#2nd_term=precipitation
        w2_rhs['g'][0]   = ivap*(upward*CC1['g'][0] - VV2['g'][0])        
        q2_rhs['g'][0]   = CC1['g'][0] - iDD*DD1['g'][0] # ivap*VV2['g'][0]#            #\nabla.{(uq2)} is added seperately 
        q1_rhs['g'][0]   = Ev['g'][0]-CC1['g'][0] #+ iDD*DD1['g'][0]    #\nabla.{(uq1)} is added seperately 
        b1_rhs['g'][0]   = - (u1['g'][0]*Db1['g'][0] + u1['g'][1]*Db1['g'][1])/a - step*(b1['g'][0]-b10['g'][0])*t_r_b1   +  step*(( CC1['g'][0]-iDD*DD1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt   +   step*RT1['g'][0]/(h1['g'][0]+H1)/tau_bu*dt #+ step*((Rad1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt                      
        b2_rhs['g'][0]   = - (u2['g'][0]*Db2['g'][0] + u2['g'][1]*Db2['g'][1])/a - step*(b2['g'][0]-b20['g'][0])*t_r_b2   +  step*((-CC1['g'][0]+iDD*DD1['g'][0]))/(h2['g'][0]+H2)/tau_bu*dt  +   step*RT2['g'][0]/(h2['g'][0]+H2)/tau_bu*dt # +   step*(Rad2['g'][0])/(h2['g'][0]+H2)/tau_bu*dt         


              
    if super_fast == 1:
        u1_rhs['g'][0]   = - (u1['g'][0]*Du1['g'][0] + u1['g'][1]*Du1['g'][2])/a - t_drag_inv*u1['g'][0] - 0.5*(h1['g'][0]*Db1['g'][0])/a-(b1['g'][0]*(Dh1['g'][0]+itopo*Dh10['g'][0]))/a - (b1['g'][0]*Dh2['g'][0])/a -  step*(u1['g'][0]-u10['g'][0])*t_r_u1
        u1_rhs['g'][1]   = - (u1['g'][0]*Du1['g'][1] + u1['g'][1]*Du1['g'][3])/a - t_drag_inv*u1['g'][1] - 0.5*(h1['g'][0]*Db1['g'][1])/a-(b1['g'][0]*Dh1['g'][1]+itopo*Dh10['g'][1])/a - (b1['g'][0]*Dh2['g'][1])/a  - step*(u1['g'][1]-u10['g'][1])*t_r_u1
        u2_rhs['g'][0]   = R2['g'][0] - (u2['g'][0]*Du2['g'][0] + u2['g'][1]*Du2['g'][2])/a - t_drag_inv*u2['g'][0] - 0.5*(h2['g'][0]*Db2['g'][0])/a-(b2['g'][0]*(Dh2['g'][0]+itopo*Dh10['g'][0]))/a - (h1['g'][0]*Db1['g'][0])/a - (b1['g'][0]*Dh1['g'][0])/a -  step*(u2['g'][0]-u20['g'][0])*t_r_u2
        u2_rhs['g'][1]   = R2['g'][1]  - (u2['g'][0]*Du2['g'][1] + u2['g'][1]*Du2['g'][3])/a - t_drag_inv*u2['g'][1] - 0.5*(h2['g'][0]*Db2['g'][1])/a-(b2['g'][0]*(Dh2['g'][1]+itopo*Dh10['g'][1]))/a - (h1['g'][0]*Db1['g'][1])/a - (b1['g'][0]*Dh1['g'][1])/a - step*(u2['g'][1]-u20['g'][1])*t_r_u2 
        h1_rhs['g'][0]   = (1.0/(b1['g'][0]+B1))  * (ep1*(-CC1['g'][0]+iDD*DD1['g'][0])+((gamma_RT1-1.0)*RT1['g'][0])) #+((gamma_NC1-1.0)*Rad1['g'][0])           
        h2_rhs['g'][0]   = (1.0/(b2['g'][0]+B2))  * (ep1*( CC1['g'][0]-iDD*DD1['g'][0])+((gamma_RT2-1.0)*RT2['g'][0]))#+((gamma_NC2-1.0)*Rad2['g'][0])
        w1_rhs['g'][0]   = ivap*((1.0-upward)*CC1['g'][0] - Prec1['g'][0])#2nd_term=precipitation
        w2_rhs['g'][0]   = ivap*(upward*CC1['g'][0] - VV2['g'][0])        
        q2_rhs['g'][0]   = ivap*VV2['g'][0]#CC1['g'][0] - iDD*DD1['g'][0]             #\nabla.{(uq2)} is added seperately 
        q1_rhs['g'][0]   = Ev['g'][0]-CC1['g'][0] #+ iDD*DD1['g'][0]    #\nabla.{(uq1)} is added seperately 
        b1_rhs['g'][0]   = - (u1['g'][0]*Db1['g'][0] + u1['g'][1]*Db1['g'][1])/a - step*(b1['g'][0]-b10['g'][0])*t_r_b1   +  step*(( CC1['g'][0]-iDD*DD1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt   +   step*RT1['g'][0]/(h1['g'][0]+H1)/tau_bu*dt #+ step*((Rad1['g'][0]))/(h1['g'][0]+H1)/tau_bu*dt                      
        b2_rhs['g'][0]   = - (u2['g'][0]*Db2['g'][0] + u2['g'][1]*Db2['g'][1])/a - step*(b2['g'][0]-b20['g'][0])*t_r_b2   +  step*((-CC1['g'][0]+iDD*DD1['g'][0]))/(h2['g'][0]+H2)/tau_bu*dt  +   step*RT2['g'][0]/(h2['g'][0]+H2)/tau_bu*dt # +   step*(Rad2['g'][0])/(h2['g'][0]+H2)/tau_bu*dt  
                
# --------------------- *** End of adding forces *** ----------------#

    '''
    # Set an upper limit for numerical stability in w1
    ww1                 = w1['g'][0]
    ws_limit            = 10.0*wcritical       #4
    ww1[ws_limit < ww1] = ws_limit
    w1['g'][0]          = ww1
    '''
       
    uh1['g'][0]      = u1['g'][0]*(h1['g'][0]) #-itopo*h10['g'][0])
    uh1['g'][1]      = u1['g'][1]*(h1['g'][0]) #-itopo*h10['g'][0])
    uh2['g'][0]      = u2['g'][0]*h2['g'][0]
    uh2['g'][1]      = u2['g'][1]*h2['g'][0]  
    #uhf1['g'][0]    = u1['g'][0]*(h1['g'][0]+H1)
    #uhf1['g'][1]    = u1['g'][1]*(h1['g'][0]+H1)
    uq1['g'][0]      = u1['g'][0]*q1['g'][0]
    uq1['g'][1]      = u1['g'][1]*q1['g'][0]
    uw1['g'][0]      = u1['g'][0]*w1['g'][0]
    uw1['g'][1]      = u1['g'][1]*w1['g'][0]
    uw2['g'][0]      = u2['g'][0]*w2['g'][0]
    uw2['g'][1]      = u2['g'][1]*w2['g'][0]     
    uq2['g'][0]      = u2['g'][0]*q2['g'][0]
    uq2['g'][1]      = u2['g'][1]*q2['g'][0]     
    #uqf2['g'][0]    = u2['g'][0]*(q2['g'][0])
    #uqf2['g'][1]    = u2['g'][1]*(q2['g'][0])   
    #uwf1['g'][0]    = u1['g'][0]*(w1['g'][0])
    #uwf1['g'][1]    = u1['g'][1]*(w1['g'][0])
    #uwf2['g'][0]    = u2['g'][0]*(w2['g'][0])
    #uwf2['g'][1]    = u2['g'][1]*(w2['g'][0])      
    #uqf1['g'][0]    = u1['g'][0]*(q1['g'][0]+Q01)
    #uqf1['g'][1]    = u1['g'][1]*(q1['g'][0]+Q01)
    # ub1['g'][0]    = u1['g'][0]*b1['g'][0]
    # ub1['g'][1]    = u1['g'][1]*b1['g'][0]
    # ubf1['g'][0]   = u1['g'][0]*(b1['g'][0]+B1)
    # ubf1['g'][1]   = u1['g'][1]*(b1['g'][0]+B1)
    
    divuh1.layout    = 'c'
    divuh2.layout    = 'c'
    divuq2.layout    = 'c'    
    divuq1.layout    = 'c'
    divuw1.layout    = 'c'
    divuw2.layout    = 'c'     
    # divub1.layout  = 'c'
    # divuhf1.layout = 'c'
    # divuqf1.layout = 'c'
    # divuqf2.layout = 'c'       
    # divubf1.layout = 'c'
    
    for m in range(m_start,m_end+1):
        md = m - m_start
        S.div(m, 1, uh1['c'][md], divuh1['c'][md])
        S.div(m, 1, uh2['c'][md], divuh2['c'][md])
        S.div(m, 1, uq1['c'][md], divuq1['c'][md])
        S.div(m, 1, uq2['c'][md], divuq2['c'][md])   
        S.div(m, 1, uw1['c'][md], divuw1['c'][md])
        S.div(m, 1, uw2['c'][md], divuw2['c'][md])               
        # S.div(m, 1, ub1['c'][md], divub1 ['c'][md])
        # S.div(m, 1, uh1['c'][md], divuhf1['c'][md])
        # S.div(m, 1, uq1['c'][md], divuqf1['c'][md])
        # S.div(m, 1, uq2['c'][md], divuqf2['c'][md])                 
        h1_rhs['c'][md] -= divuh1['c'][md]/a
        h2_rhs['c'][md] -= divuh2['c'][md]/a
        q2_rhs['c'][md] -= divuq2['c'][md]/a        
        q1_rhs['c'][md] -= divuq1['c'][md]/a
        w1_rhs['c'][md] -= divuw1['c'][md]/a
        w2_rhs['c'][md] -= divuw2['c'][md]/a

    RHS.pack(u1_rhs, h1_rhs, u2_rhs, h2_rhs, q1_rhs, b1_rhs, b2_rhs, q2_rhs, w1_rhs, w2_rhs)
    state_vector.pack(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)
# end of nonlinear function   

# Output variables are nondimensional; the corresponding dimensional scales are provided below. Dimensional scales may vary depending on the
# nondimensionalization method.
# The dimensional scale of velocity on the equatorial beta plane is beta * L_d^2, where L_d represents the equatorial Rossby deformation radius.
# Layer 1 corresponds to the lower layer, and layer 2 to the upper layer.
# The variable q represents the bulk of humidity, with a scale of [L_v * g * H1 / (C_p * theta_s * 45)]. However, it can be rescaled to the vertically averaged
# specific humidity with a scale of [(1/45) Kg/Kg].
# The time scale T is given by T = 1 / (beta * L_d^2) if x and y are nondimensionalized by the equatorial Rossby deformation radius.
# RT stands for Thermal forcing due to Radiative Transfer.



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
#ncout.AddVariable('Press1',    'Pressure layer 1',                            '[Pa]',                                     miss_val=None, axes=axes_tyx)
#ncout.AddVariable('Press2',    'Pressure layer 2',                            '[Pa]',                                     miss_val=None, axes=axes_tyx)

# ncout.AddVariable('R2th',     'long name of R2th',                           'unknown units',                            miss_val=None, axes=axes_tyx)
# ncout.AddVariable('R2ph',     'long name of R2ph',                           'unknown units',                            miss_val=None, axes=axes_tyx)
# ncout.AddVariable('divuhf1',  'long name of divuhf1',                        'unknown units',                            miss_val=None, axes=axes_tyx)
# ncout.AddVariable('divuqf1',  'long name of divuqf1',                        'unknown units',                            miss_val=None, axes=axes_tyx)
# ncout.AddVariable('divuwf1',  'long name of divuwf1',                        'unknown units',                            miss_val=None, axes=axes_tyx)
# ncout.AddVariable('divuwf2',  'long name of divuwf2',                        'unknown units',                            miss_val=None, axes=axes_tyx)
# ncout.AddVariable('divubf1',  'long name of divubf1',                        'unknown units',                            miss_val=None, axes=axes_tyx)
ncout.AddVariable('menthalpy',  'menthalpy = h1+H1-ep1*q1-ep1*Q01',            'unknown units',                            miss_val=None, axes=axes_tyx)
ncout.AddVariable('insolation', 'daily mean insolation',                       'W/m2?',                                    miss_val=None, axes=axes_ty) 
ncout.WriteVar('h10',    h10['g'][0],    t*T_scale)
ncout.WriteVar('ocean',   Oc['g'][0],    t*T_scale) 
np.seterr(all='raise')
mainloop_time = time.time()


print0()
print0()
print0('Informations:')
print0()
print0('Variables            Description')
print0('-------------------------------------------------------')
print0('h1,h2                Layers Pseudo-height')
print0('q1                   Bulk of Specific humidity')
print0('m_enthalpy1          Enthalpy')
print0('b1,b2                Buoyancy')
print0('u1,u2                Meridional velocity')
print0('w1                   Bulk of Precipitable Water')
print0('E1                   Sea surface evaoporation')
print0('C1                   CLWC (Cloud Liquid Water Content)')
print0('D1                   Downdraft to layer 1 (balanced)')
print0('nd                   Non dimensional')
print0('-------------------------------------------------------')

print0()
print0()
print0('Simulation progress:')

ccc = 0

# Main loop of actual sequence to iterate over
for i in range(n_iterations):

    if i % n_output == 0:

        state_vector.unpack(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)

        # for m in range(m_start,m_end+1):
            # md = m - m_start
            # (start_index,end_index,spin) = S.tensor_index(m,1)
            # om['c'][md] = 1j*(S.op('k-',m,1).dot(u['c'][md][start_index[0]:end_index[0]]) - S.op('k+',m,-1).dot(u['c'][md][start_index[1]:end_index[1]]))

        # Write NetCDFOutput in all ranks at output time.
        # print("u1th ", u1['g'][0])
        ncout.WriteVar('u1th',           u1['g'][0], t*T_scale)
        ncout.WriteVar('u1ph',           u1['g'][1], t*T_scale)
        ncout.WriteVar('h1',             h1['g'][0], t*T_scale)
        ncout.WriteVar('q1',             q1['g'][0], t*T_scale)
        ncout.WriteVar('q2',             q2['g'][0], t*T_scale)
        ncout.WriteVar('w1',             w1['g'][0], t*T_scale)
        ncout.WriteVar('w2',             w2['g'][0], t*T_scale)
        ncout.WriteVar('b1',             b1['g'][0], t*T_scale)
        ncout.WriteVar('b2',             b2['g'][0], t*T_scale)
        ncout.WriteVar('CC1',           CC1['g'][0], t*T_scale)
        ncout.WriteVar('Prec1',       Prec1['g'][0], t*T_scale)
        ncout.WriteVar('DD1',           DD1['g'][0], t*T_scale)
        ncout.WriteVar('DDr',           DDr['g'][0], t*T_scale)
        ncout.WriteVar('Ev',             Ev['g'][0], t*T_scale)
        ncout.WriteVar('u2th',           u2['g'][0], t*T_scale)
        ncout.WriteVar('u2ph',           u2['g'][1], t*T_scale)
        ncout.WriteVar('h2',             h2['g'][0], t*T_scale)
        ncout.WriteVar('RT1',           RT1['g'][0], t*T_scale)
        ncout.WriteVar('RT2',           RT2['g'][0], t*T_scale)        
        #ncout.WriteVar('Press1',    Press1['g'][0], t*T_scale)     
        #ncout.WriteVar('Press2',    Press2['g'][0], t*T_scale)            
        # ncout.WriteVar('R2th',         R2['g'][0], t*T_scale)
        # ncout.WriteVar('R2ph',         R2['g'][1], t*T_scale)
        # ncout.WriteVar('divuhf1', divuhf1['g'][0], t*T_scale)
        # ncout.WriteVar('divuqf1', divuqf1['g'][0], t*T_scale)
        # ncout.WriteVar('divuwf1', divuwf1['g'][0], t*T_scale)
        # ncout.WriteVar('divuwf2', divuwf2['g'][0], t*T_scale)
        # ncout.WriteVar('divubf1', divubf1['g'][0], t*T_scale)
        
        ment_local       = h1['g'][0] + H1 - ep1*q1['g'][0] - ep1*Q01
        ncout.WriteVar('menthalpy',      ment_local, t*T_scale)

        # gather full data to output
        u1th_global      = comm.gather(u1['g'][0],      root=0)
        u1ph_global      = comm.gather(u1['g'][1],      root=0)
        h1_global        = comm.gather(h1['g'][0],      root=0)
        q1_global        = comm.gather(q1['g'][0],      root=0)
        q2_global        = comm.gather(q2['g'][0],      root=0)   
        w1_global        = comm.gather(w1['g'][0],      root=0)
        w2_global        = comm.gather(w2['g'][0],      root=0)                       
        b1_global        = comm.gather(b1['g'][0],      root=0)
        b2_global        = comm.gather(b2['g'][0],      root=0)        
        CC1_global       = comm.gather(CC1['g'][0],     root=0)
        VV1_global       = comm.gather(VV1['g'][0],     root=0)        
        VV2_global       = comm.gather(VV2['g'][0],     root=0) 
        Prec1_global     = comm.gather(Prec1['g'][0],   root=0)        
        DD1_global       = comm.gather(DD1['g'][0],     root=0)        
        DDr_global       = comm.gather(DDr['g'][0],     root=0)                
        Ev_global        = comm.gather(Ev['g'][0],      root=0)
        u2th_global      = comm.gather(u2['g'][0],      root=0)
        u2ph_global      = comm.gather(u2['g'][1],      root=0)
        h2_global        = comm.gather(h2['g'][0],      root=0)
        RT1_global       = comm.gather(RT1['g'][0],     root=0)  
        RT2_global       = comm.gather(RT2['g'][0],     root=0)    
        #Press1_global   = comm.gather(Press1['g'][0],  root=0)
        #Press2_global   = comm.gather(Press2['g'][0],  root=0)                         
        # R2th_global    = comm.gather(R2['g'][0],      root=0)
        # R2ph_global    = comm.gather(R2['g'][1],      root=0)
        # divuhf1_global = comm.gather(divuhf1['g'][0], root=0)
        # divuqf1_global = comm.gather(divuqf1['g'][0], root=0)
        # divuqf2_global = comm.gather(divuqf2['g'][0], root=0)              
        # divubf1_global = comm.gather(divubf1['g'][0], root=0)
        
                        
        if rank == 0:
            # Save data
            u1ph_global     = np.hstack(u1ph_global)
            u1th_global     = np.hstack(u1th_global)
            h1_global       = np.hstack(h1_global)
            q1_global       = np.hstack(q1_global)      
            b1_global       = np.hstack(b1_global)
            b2_global       = np.hstack(b2_global)            
            u2ph_global     = np.hstack(u2ph_global)
            u2th_global     = np.hstack(u2th_global)
            h2_global       = np.hstack(h2_global)
            RT1_global      = np.hstack(RT1_global)     
            RT2_global      = np.hstack(RT2_global)                               
            # R2ph_global   = np.hstack(R2ph_global)
            # R2th_global   = np.hstack(R2th_global)
            CC1_global      = np.hstack(CC1_global)
            # VV1_global    = np.hstack(VV1_global)   
            # VV2_global    = np.hstack(VV2_global)  
            Prec1_global    = np.hstack(Prec1_global)            
            DD1_global      = np.hstack(DD1_global)    
            DDr_global      = np.hstack(DDr_global)                                
            Ev_global       = np.hstack(Ev_global)
            q2_global       = np.hstack(q2_global)
            w1_global       = np.hstack(w1_global)
            w2_global       = np.hstack(w2_global)                  
            menthalpy       = h1_global+H1-ep1*q1_global-ep1*Q01
            
            #Press1_global  = np.hstack(Press1_global)
            #Press2_global  = np.hstack(Press2_global)
            
            output_filename = '{}/{}{:05d}'.format(output_folder, output_file, file_num)  # example: ('output', 'output_', 1)
            
            np.savez_compressed(output_filename,
                     # p       = p_global, 
                     # om      = om_global, 
                     # vph     = vph_global, 
                     # vth     = vth_global,
                     h1        = h1_global, 
                     u1ph      = u1ph_global, 
                     u1th      = u1th_global,
                     h2        = h2_global, 
                     RT1       = RT1_global,
                     RT2       = RT2_global,                    
                     u2ph      = u2ph_global, 
                     u2th      = u2th_global,
                     q1        = q1_global, 
                     # R2ph    = R2ph_global, 
                     # R2th    = R2th_global,
                     # divubf1 = divubf1_global, 
                     b1        = b1_global, 
                     b2        = b2_global, 
                     q2        = q2_global, 
                     w1        = w1_global, 
                     w2        = w2_global,
                     # divuqf2 = divuqf2_global,
                     menthalpy = menthalpy,
                     # divuhf1 = divuhf1_global, 
                     # divuqf1 = divuqf1_global,
                     CC1       = CC1_global, 
                     # VV1     = VV1_global, 
                     # VV2     = VV2_global,
                     Prec1     = Prec1_global,
                     DD1       = DD1_global,
                     DDr       = DDr_global, 
                     Ev        = Ev_global, 
                     #Press1   = Press1_global, 
                     #Press2   = Press2_global, 
                     t         = np.array([t]), 
                     lamda     = lamda[:,0], 
                     theta     = theta_global)
            file_num += 1

            # Print iteration and maximum vorticity
            
            # step 
            step_h1_min      = np.min(h1_global)
            step_h1_max      = np.max(h1_global)
            step_h2_min      = np.min(h2_global)
            step_h2_max      = np.max(h2_global)
            step_q1_max      = np.max(q1_global)
            step_m_enthalpy1 = np.min(menthalpy)
            step_b1_max      = np.max(b1_global)
            step_b2_max      = np.max(b2_global)
            step_u1_max      = np.max(u1ph_global)
            step_w1_max      = np.max(w1_global)
            #step_p1_max     = np.max(Press1_global)
            #step_p2_max     = np.max(Press2_global)
            #step_p1_min     = np.min(Press1_global)
            #step_p2_min     = np.min(Press2_global)
            # totals
            total_q1         = q1_global.sum()
            total_E1         = Ev_global.sum()
            total_C1         = CC1_global.sum()
            total_D1         = DD1_global.sum()
            # total_Press1   = Press1_global.sum()
            # total_Press2   = Press2_global.sum()
                                    
            text_print_1  = '--> Iteration =  {0:7} | Time {1:7} [days] = {2:7} [nondim. days]'.format(i, np.round(t*T_scale,4), np.round(t,4))
            
            text_print_2  = '--> Step:  h1_min {0:6.3f} | h1_max {1:6.3f} | h2_min {2:6.3f} | h2_max {3:6.3f} | q1_max {4:6.3f} | m_enthalpy1 {5:5.3f} | b1_max {6:5.3f} | b2_max {7:5.3f} | u1_max {8:5.3f} | w1_max {9:5.3f} ||| Totals: q1 {10:6.1f} | E1 {11:5.3f} | C1 {12:5.3f} | D1 {13:5.3f}'.\
                                format(step_h1_min, step_h1_max, step_h2_min, step_h2_max, step_q1_max, step_m_enthalpy1, step_b1_max, step_b2_max, step_u1_max, step_w1_max, total_q1, total_E1, total_C1, total_D1)
                               
            text_print_3  = '--> Output = {}.npz ({} | running for {})'.format(output_filename, time.asctime(), str(timedelta(seconds=time.time()-start_time))[:-7]) 
                        
            # --> Iteration =      285 | Time  1.2435 [days] =  1.4445 [nondim. days]
            # --> Step:  h1_min -0.056 | h1_max  0.065 | h2_min -0.074 | h2_max  0.068 | q1_max  0.061 | m_enthalpy1 0.544 | b1_max 0.122 | b2_max 0.092 | u1_max 0.068 | w1_max 0.000 ||| Totals: q1   73.5 | E1 0.000 | C1 0.000 | D1 0.000
            # --> Output = output/output_00006.npz (Wed Jun 12 17:50:54 2024 | running for 0:03:07)
            
            print()
            print(text_print_1)
            print(text_print_2)
            print(text_print_3)
                       
    
    ##########################################################################
    ############################ BUMP FORCING ################################
    ##########################################################################

    if external_forcing != 0 and int(t*T_scale) == 5 and ccc == 0:  # Start the Bump forcing at N real days
        
        # defined at the beginning --> epsilon = external_forcing_epsilon  # from config.py
       
        b1_bump = epsilon * H1* np.sqrt(2*np.exp(1)) * 2**(1./spara)/spara * psi    # Bump perturbation 
        b1_bump = np.where(b1_bump <= -0.01, 0, b1_bump) 

        b2_bump = epsilon * H2* np.sqrt(2*np.exp(1)) * 2**(1./spara)/spara * psi    # Bump perturbation 
        b2_bump = np.where(b2_bump <= -0.01, 0, b2_bump) 
        
        if external_forcing == 'baroclinic': 
            b1['g'][0]  += b1_bump      # adding perturbation to lower layer
        if external_forcing == 'barotropic': 
            b1['g'][0]  += b1_bump      # adding perturbation to lower layer
            b2['g'][0]  += b2_bump       
        
        state_vector.pack(u1, h1, u2, h2, q1, b1, b2, q2, w1, w2)
        
        ccc += 1
        
        print0()
        print0('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print0('+   BUMP ADDED ==>   epsilon {:5} sh_E {:6} sh_N {:6}+'.format(str(epsilon), str(sh_E), str(sh_N)))
        print0('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        
    ##########################################################################
    ##########################################################################
    ##########################################################################
      
            
    nonlinear(state_vector, RHS)
    timestepper.step(dt, state_vector, S, L, M, P, RHS, LU)
    t       += dt


#    # imposing that the m=0 mode of u,h,c are purely real
#    if i % 100 == 1:
#        state_vector.unpack(u,h,c)
#        u.require_grid_space()
#        u.require_coeff_space()
#        h.require_grid_space()
#        h.require_coeff_space()
#        c.require_grid_space()
#        c.require_coeff_space()
#        state_vector.pack(u,h,c)


ncout.close()
end_time = time.time()
print0()
print0('simulated:          ', i, 'timesteps -->', (t*T_scale), 'days =', t, 'nd' )
print0('total time:         ', end_time-start_time)
print0('init time:          ', mainloop_time-start_time)
print0('mainloop time:      ', end_time-mainloop_time)



