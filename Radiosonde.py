#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 14:49:17 2020

@author: ncp532
"""

# Drawing packages
import matplotlib.pyplot as plt
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

# Data handing packages
import numpy as np
import pandas as pd
from scipy import signal, stats
import netCDF4 as nc
import xarray as xr
from netCDF4 import Dataset

# Date and Time handling package
import datetime as dt
from datetime import datetime,time, timedelta		# functions to handle date and time
from matplotlib.colors import BoundaryNorm

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# Radiosonde
RS_2017 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/Radiosonde/ARM_Radiosonde_20171101_20180324.csv',   index_col=0)

#------------------------------------------------------------------------------
# SPLIT THE RADIOSONDE ARRY INTO INDIVIDUAL SONDES

# Set the date
RS_2017.index   = (pd.to_datetime(RS_2017.index,   dayfirst=True))

#------------------------------
# Variables
#------------------------------
Time  = RS_2017.index        # UT Time (YYYY-MM-DD HH:MM:SS)
Alt   = RS_2017['alt']       # Altitude (m)
Pres  = RS_2017['pres']      # Pressure (hPa)
Tdry  = RS_2017['tdry']      # Dry bulb temperature (C)
Wspd  = RS_2017['wspd']      # Wind speed (m/s)
Wdir  = RS_2017['deg']       # Wind direction (deg)
RH    = RS_2017['rh']        # Relative humidity (%)
Uwind = RS_2017['u_wind']    # Eastward wind component (m/s)
Vwind = RS_2017['v_wind']    # Northward wind component (m/s)

#------------------------------
# Seperate the individual sondes
#------------------------------
 # array of sub-arrays (starts with first value)
UTTime    = [[Time[0]]]
Altitude  = [[Alt[0]]]
Pressure  = [[Pres[0]]]
TempDry   = [[Tdry[0]]]
WindSpeed = [[Wspd[0]]]
WindDir   = [[Wdir[0]]]
RelHum    = [[RH[0]]]
EastWind  = [[Uwind[0]]]
NorthWind = [[Vwind[0]]]

# go through each element based on the altitude
for i in range(1, len(Alt)):
    # If altitude is larger than previous altitude
    if Alt[i - 1] <= Alt[i]:
        # Add the variables to the last sub-array
        Altitude[len(Altitude)   - 1].append(Alt[i])
        Pressure[len(Pressure)   - 1].append(Pres[i])
        UTTime[len(UTTime)       - 1].append(Time[i])
        TempDry[len(TempDry)     - 1].append(Tdry[i])
        WindSpeed[len(WindSpeed) - 1].append(Wspd[i])
        WindDir[len(WindDir)     - 1].append(Wdir[i])
        RelHum[len(RelHum)       - 1].append(RH[i])
        EastWind[len(EastWind)   - 1].append(Uwind[i])
        NorthWind[len(NorthWind) - 1].append(Vwind[i])
    # If the altitude isn't larger than previous altitude
    else:
        # Add the variables to a new sub-array
        Altitude.append([Alt[i]])
        Pressure.append([Pres[i]])
        UTTime.append([Time[i]])
        TempDry.append([Tdry[i]])
        WindSpeed.append([Wspd[i]])
        WindDir.append([Wdir[i]])
        RelHum.append([RH[i]])
        EastWind.append([Uwind[i]])
        NorthWind.append([Vwind[i]])

# Convert the results to Pandas DataFrames 
Altitude  = pd.DataFrame(Altitude)
Pressure  = pd.DataFrame(Pressure)
UTTime    = pd.DataFrame(UTTime)
TempDry   = pd.DataFrame(TempDry)
WindSpeed = pd.DataFrame(WindSpeed)
WindDir   = pd.DataFrame(WindDir)
RelHum    = pd.DataFrame(RelHum)
EastWind  = pd.DataFrame(EastWind)
NorthWind = pd.DataFrame(NorthWind)

# Convert UTTime to an index for the other DataFrames
IndexTime = UTTime.T.loc[0]

# Levels
Levels = np.array(UTTime.columns)

# Set the DataFrame index to IndexTime
UTTime.index    = IndexTime
Altitude.index  = IndexTime
Pressure.index  = IndexTime
TempDry.index   = IndexTime
WindSpeed.index = IndexTime
WindDir.index   = IndexTime
RelHum.index    = IndexTime
EastWind.index  = IndexTime
NorthWind.index = IndexTime

# Resample the dataset to 1-hour resolution
Altitude  = Altitude.resample('60T').mean()
Pressure  = Pressure.resample('60T').mean()
TempDry   = TempDry.resample('60T').mean()
WindSpeed = WindSpeed.resample('60T').mean()
WindDir   = WindDir.resample('60T').mean()
RelHum    = RelHum.resample('60T').mean()
EastWind  = EastWind.resample('60T').mean()
NorthWind = NorthWind.resample('60T').mean()

# Interpolate the values within the dataset
Altitude  = Altitude.interpolate(method='time')
Pressure  = Pressure.interpolate(method='time')
TempDry   = TempDry.interpolate(method='time')
WindSpeed = WindSpeed.interpolate(method='time')
WindDir   = WindDir.interpolate(method='time')
RelHum    = RelHum.interpolate(method='time')
EastWind  = EastWind.interpolate(method='time')
NorthWind = NorthWind.interpolate(method='time')
IndexTime = Altitude.index

#------------------------------
# Save as a NetCDF File
#------------------------------

#Define a filename (fn)
fn = '/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/Radiosonde/test.nc'
# Call dataset (ds) and sepecify write mode (w)
ds = nc.Dataset(fn, 'w', format='NETCDF4')

# Add dimensions
DateTime = ds.createDimension('DateTime', None)
Level    = ds.createDimension('Level',    None)

# Add NetCDF variables
times     = ds.createVariable('DateTime',  'f4', ('DateTime',))
levels    = ds.createVariable('Level',     'f4', ('Level',))
altitude  = ds.createVariable('altitude',  'f4', ('DateTime', 'Level',))
pressure  = ds.createVariable('pressure',  'f4', ('DateTime', 'Level',))
tempdry   = ds.createVariable('tempdry',   'f4', ('DateTime', 'Level',))
windspeed = ds.createVariable('windspeed', 'f4', ('DateTime', 'Level',))
winddir   = ds.createVariable('winddir',   'f4', ('DateTime', 'Level',))
relhum    = ds.createVariable('relhum',    'f4', ('DateTime', 'Level',))
eastwind  = ds.createVariable('eastwind',  'f4', ('DateTime', 'Level',))
northwind = ds.createVariable('northwind', 'f4', ('DateTime', 'Level',))

# Define the units for the NetCDF variables
times.units     = 'YYYY-MM-DD HH:MM:SS'
levels.units    = 'Unknown'
altitude.units  = 'm'
pressure.units  = 'hPa'
tempdry.units   = 'C'
windspeed.units = 'm/s'
winddir.units   = 'deg'
relhum.units    = '%'
eastwind.units  = 'm/s'
northwind.units = 'm/s'

# Assign time and level values
times[:]  = IndexTime
levels[:] = Levels

# Assign data values to the variables
altitude.units  = Altitude
pressure.units  = Pressure
tempdry.units   = TempDry
windspeed.units = WindSpeed
winddir.units   = WindDir
relhum.units    = RelHum
eastwind.units  = EastWind
northwind.units = NorthWind

# Close the dataset
ds.close()
