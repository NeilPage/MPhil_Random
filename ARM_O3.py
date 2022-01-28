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

ARM = Dataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/Original_O3/maraoso3M1.b1.20180324.000000.nc')
ARM2 = MFDataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/maraoscoM1.a1.201*.nc')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES
Time = ARM.variables[u'time'][:] # Time (seconds since midnight 2017/10/18)
O3 = ARM.variables[u'o3'][:] # Ozone (O3)

#------------------------------------------------------------------------------
# SET DATETIME
d0=datetime(2018,3,24,0,0,0) # this date corresponds to the date of the arm file
date = []
for t in Time:
    sec=timedelta(seconds=int(t))
    date.append(d0+sec)

#------------------------------------------------------------------------------
# BUILD A DATAFRAME AND SAVE AS A .CSV
dfO3 = np.column_stack((date, O3))
dfO3 = pd.DataFrame.from_dict(dfO3)
dfO3.columns = ['Date', 'O3_(ppb)']
dfO3.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/test.csv')
