#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:42:11 2019

@author: ncp532
"""

# FILE SYSTEM PACKAGES
from netCDF4 import Dataset			# function used to open single netcdf file
from netCDF4 import MFDataset		# function used to open multiple netcdf files
import xarray as xr                 # xarray to open multiple netcdf files    

# DATE AND TIME HANDLING PACKAGES
from datetime import datetime,timedelta		# functions to handle date and time

# DATA HANDLING PACKAGES
import numpy as np					    # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.
import pandas as pd

#------------------------------------------------------------------------------
# DEFINE THE DATASET

ARM = Dataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/Original_CO_N2O/maraoscoM1.a1.20171018.234750.nc')
#ARM2 = MFDataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/maraoscoM1.a1.201*.nc')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES
Time = ARM.variables[u'time'][:] # Time (seconds since midnight 2017/10/18)
CO = ARM.variables[u'co'][:] # Carbon Monoxide (CO)
N2O = ARM.variables[u'n2o'][:] # Nitrous Oxide (N2O)
CO_Dry = ARM.variables[u'co_dry'][:] # Carbon Monoxide Dry (CO)
N2O_Dry = ARM.variables[u'n2o_dry'][:] # Nitrous Oxide Dry (N2O)

#------------------------------------------------------------------------------
# SET DATETIME
d0=datetime(2017,10,18,23,47,50) # this date corresponds to the date of the arm file
date = []
for t in Time:
    sec=timedelta(seconds=int(t))
    date.append(d0+sec)

#------------------------------------------------------------------------------
# BUILD A DATAFRAME AND SAVE AS A .CSV
dfARM = np.column_stack((date, CO, N2O, CO_Dry, N2O_Dry))
dfARM = pd.DataFrame.from_dict(dfARM)
dfARM.columns = ['Date', 'CO','N2O','CO_Dry','N2O_Dry']
dfARM.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM_20171018.csv')
