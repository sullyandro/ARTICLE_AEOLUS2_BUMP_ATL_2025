"""
aeolus2test.py

A mininmal stand-alone main program to test the Aeolus2 functionality
that is implemented in aeolus2.py
"""

import sys
import os

import numpy as np

from mpi4py import MPI

# aeolus2.py implements the real work.
from aeolus2 import \
    aeolus2_init_grid, aeolus2_init, \
    aeolus2_restart, aeolus2_finish, \
    aeolus2_get_stock_pe, \
    aeolus2_get_bottom_mass, aeolus2_get_bottom_wind, \
    aeolus2_update_down, aeolus2_update_up


# start time of experiment
Time_init_year = 1
Time_init_month = 1
Time_init_day = 1
Time_init_seconds = 0
# start time of this run, within the experiment
Time_days = 1
Time_seconds = 0

Time_step_days = 0
Time_step_seconds = 3600

without_topo = False
update_land_frac_always = False

nr_tracers = 1    # number of tracers to exchange with coupler
q_ind = 1         # index of humidity tracer
co2_ind = -1      # index of CO2 tracer


ics, ice, jcs, jce, lons, lonb, lats_global, latb_global = \
   aeolus2_init_grid(nlons_global=192, nlats_global=96, lon_0_360=False, comm_in=MPI.COMM_WORLD)

ni = ice-ics
nj = jce-jcs

HORO = np.zeros(ni*nj, dtype=np.float64)
SIGORO = np.zeros(ni*nj, dtype=np.float64)
land_frac = np.zeros(ni*nj, dtype=np.float64)
area = np.zeros(ni*nj, dtype=np.float64)

HORO = np.reshape(HORO,(ni,nj),order='F')
SIGORO = np.reshape(SIGORO,(ni,nj),order='F')
land_frac = np.reshape(land_frac,(ni,nj),order='F')
area = np.reshape(area,(ni,nj),order='F')

aeolus2_init((Time_init_year, Time_init_month, Time_init_day, Time_init_seconds,
              Time_days, Time_seconds,
              Time_step_days, Time_step_seconds,
              without_topo, update_land_frac_always,
              nr_tracers, q_ind, co2_ind,
              HORO, SIGORO, land_frac, area))


aeolus2_finish(Time_seconds, Time_days)

MPI.Finalize()
