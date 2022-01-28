#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:51:48 2020

@author: ncp532
"""

# FILE SYSTEM PACKAGES
from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
import xarray as xr

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# DATA HANDLING PACKAGES
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from scipy import stats
import math

# DRAWING PACKAGES
import cartopy.crs as ccrs
from matplotlib import cm                   # imports the colormap function from matplotlib
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt             
import matplotlib.dates as mdates            
from matplotlib.ticker import MaxNLocator
import cmocean
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#------------------------------------------------------------------------------
# OBSERVATIONAL DATA

#--------------
# BrO
#--------------
BrO_V1_17   = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V1_17/all_BrO/V1_17_BrO_retrieval.csv', index_col=0)       # BrO V1 (2017/18)
BrO_V2_17   = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V2_17/all_BrO/V2_17_BrO_retrieval.csv', index_col=0)       # BrO V2 (2017/18)
BrO_V3_17   = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V3_17/all_BrO/V3_17_BrO_retrieval.csv', index_col=0)       # BrO V3 (2017/18)

BrO_V1_18   = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V1_18/all_BrO/V1_18_BrO_retrieval.csv', index_col=0)       # BrO V1 (2018/19)
BrO_V2_18   = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V2_18/all_BrO/V2_18_BrO_retrieval.csv', index_col=0)       # BrO V2 (2018/19)
BrO_V3_18   = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V3_18/all_BrO/V3_18_BrO_retrieval.csv', index_col=0)       # BrO V3 (2018/19)

BrO_SIPEXII = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/SIPEXII/all_BrO/SIPEXII_BrO_retrieval.csv', index_col=0) # BrO SIPEXII (2012)

#--------------
# Hg0
#--------------
Hg0_V1_17   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V1_Hg0_QAQC_17-18.csv', index_col=0)
Hg0_V2_17   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V2_Hg0_QAQC_17-18.csv', index_col=0)
Hg0_V3_17M  = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V3_Hg0_QAQC_17-18.csv', index_col=0)
Hg0_V3_17D  = Hg0_V3_17M

Hg0_V1_18   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/CAMMPCAN_V1_Hg0_QAQC_18-19.csv', index_col=0)
Hg0_V2_18   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/CAMMPCAN_V2_Hg0_QAQC_18-19.csv', index_col=0)
Hg0_V3_18M  = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/CAMMPCAN_V3_Hg0_QAQC_18-19.csv', index_col=0)
Hg0_V3_18D  = Hg0_V3_18M

Hg0_SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/SIPEXII_Hg_Air/SIPEXII_Hg0_QAQC_2012.csv', index_col=0)

#-------------
# Sea Ice
#-------------

SeaIce_V1_17   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2017-18/oisst-avhrr-v02r01.20171114.nc', decode_cf=False, engine='netcdf4')
SeaIce_V2_17   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2017-18/oisst-avhrr-v02r01.20171221.nc', decode_cf=False, engine='netcdf4')
SeaIce_V3_17M  = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2017-18/oisst-avhrr-v02r01.20180201.nc', decode_cf=False, engine='netcdf4')
SeaIce_V3_17D  = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2017-18/oisst-avhrr-v02r01.20180127.nc', decode_cf=False, engine='netcdf4')

SeaIce_V1_18   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20181107.nc', decode_cf=False, engine='netcdf4')
SeaIce_V2_18   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20181215.nc', decode_cf=False, engine='netcdf4')
SeaIce_V3_18M  = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20190130.nc', decode_cf=False, engine='netcdf4')
SeaIce_V3_18D  = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20190126.nc', decode_cf=False, engine='netcdf4')

SeaIce_SIPEXII = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2012/oisst-avhrr-v02r01.20120923.nc', decode_cf=False, engine='netcdf4')
SeaIce_PCAN    = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20170126_median5day.nc', decode_cf=False, engine='netcdf4')

#------------------------------------------------------------------------------
# FILTER THE BrO DATA FOR RELATIVE ERROR 

#--------------
# BrO
#--------------
# Calculate the Relative Error (>=0.6)
Filter1_BrO = BrO_V1_17['err_surf_vmr'] / BrO_V1_17['surf_vmr(ppmv)']
Filter2_BrO = BrO_V2_17['err_surf_vmr'] / BrO_V2_17['surf_vmr(ppmv)']
Filter3_BrO = BrO_V3_17['err_surf_vmr'] / BrO_V3_17['surf_vmr(ppmv)']

Filter4_BrO = BrO_V1_18['err_surf_vmr'] / BrO_V1_18['surf_vmr(ppmv)']
Filter5_BrO = BrO_V2_18['err_surf_vmr'] / BrO_V2_18['surf_vmr(ppmv)']
Filter6_BrO = BrO_V3_18['err_surf_vmr'] / BrO_V3_18['surf_vmr(ppmv)']

Filter7_BrO = BrO_SIPEXII['err_surf_vmr'] / BrO_SIPEXII['surf_vmr(ppmv)']

# Apply the filter
V1_17F      = Filter1_BrO < 0.6
BrO_V1_17   = BrO_V1_17[V1_17F]

V2_17F      = Filter2_BrO < 0.6
BrO_V2_17   = BrO_V2_17[V2_17F]

V3_17F      = Filter3_BrO < 0.6
BrO_V3_17   = BrO_V3_17[V3_17F]

V1_18F      = Filter4_BrO < 0.6
BrO_V1_18   = BrO_V1_18[V1_18F]

V2_18F      = Filter5_BrO < 0.6
BrO_V2_18   = BrO_V2_18[V2_18F]

V3_18F      = Filter6_BrO < 0.6
BrO_V3_18   = BrO_V3_18[V3_18F]

SIPEXIIF    = Filter7_BrO < 0.6
BrO_SIPEXII = BrO_SIPEXII[SIPEXIIF]

# Split V3 into Mawson and Davis
BrO_V3_17M = BrO_V3_17
BrO_V3_17D = BrO_V3_17

BrO_V3_18M = BrO_V3_18
BrO_V3_18D = BrO_V3_18

#------------------------------------------------------------------------------
# Set the date

#--------------
# BrO
#--------------
BrO_V1_17.index   = (pd.to_datetime(BrO_V1_17.index,   dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7
BrO_V2_17.index   = (pd.to_datetime(BrO_V2_17.index,   dayfirst=True) + timedelta(hours=8)) # Casey timezone is UT+8
BrO_V3_17M.index  = (pd.to_datetime(BrO_V3_17M.index,  dayfirst=True) + timedelta(hours=5)) # Mawson timezone is UT+5
BrO_V3_17D.index  = (pd.to_datetime(BrO_V3_17D.index,  dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7

BrO_V1_18.index   = (pd.to_datetime(BrO_V1_18.index,   dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7
BrO_V2_18.index   = (pd.to_datetime(BrO_V2_18.index,   dayfirst=True) + timedelta(hours=8)) # Casey timezone is UT+8
BrO_V3_18M.index  = (pd.to_datetime(BrO_V3_18M.index,  dayfirst=True) + timedelta(hours=5)) # Mawson timezone is UT+5
BrO_V3_18D.index  = (pd.to_datetime(BrO_V3_18D.index,  dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7

BrO_SIPEXII.index = (pd.to_datetime(BrO_SIPEXII.index, dayfirst=True) + timedelta(hours=8)) # SIPEXII timezone is UT+8

#--------------
# Hg0
#--------------
Hg0_V1_17.index   = pd.to_datetime(Hg0_V1_17.index,   dayfirst=True)
Hg0_V2_17.index   = pd.to_datetime(Hg0_V2_17.index,   dayfirst=True)
Hg0_V3_17M.index  = pd.to_datetime(Hg0_V3_17M.index,  dayfirst=True)
Hg0_V3_17D.index  = (pd.to_datetime(Hg0_V3_17D.index, dayfirst=True) - timedelta(hours=1))

Hg0_V1_18.index   = pd.to_datetime(Hg0_V1_18.index,   dayfirst=True)
Hg0_V2_18.index   = pd.to_datetime(Hg0_V2_18.index,   dayfirst=True)
Hg0_V3_18M.index  = pd.to_datetime(Hg0_V3_18M.index,  dayfirst=True)
Hg0_V3_18D.index  = (pd.to_datetime(Hg0_V3_18D.index, dayfirst=True) - timedelta(hours=1))

Hg0_SIPEXII.index = pd.to_datetime(Hg0_SIPEXII.index, dayfirst=True)

#------------------------------------------------------------------------------
# PASSIVATION ISSUE WITH CELL A ON VOYAGES V3_18M, V3_18D and V4_18 (FILTER DATA)

Filter1    = Hg0_V3_18M['Cart'] == "B"
Hg0_V3_18M = Hg0_V3_18M[Filter1]

Filter2    = Hg0_V3_18D['Cart'] == "B"
Hg0_V3_18D = Hg0_V3_18D[Filter2]

#------------------------------------------------------------------------------
# RESAMPLE THE Hg0 DATASETS TO 1-HOUR TIME RESOLUTION

#--------------
# Hg0
#--------------
Hg0_V1_17   = Hg0_V1_17.resample('60T').mean()
Hg0_V2_17   = Hg0_V2_17.resample('60T').mean()
Hg0_V3_17M  = Hg0_V3_17M.resample('60T').mean()
Hg0_V3_17D  = Hg0_V3_17D.resample('60T').mean()

Hg0_V1_18   = Hg0_V1_18.resample('60T').mean()
Hg0_V2_18   = Hg0_V2_18.resample('60T').mean()
Hg0_V3_18M  = Hg0_V3_18M.resample('60T').mean()
Hg0_V3_18D  = Hg0_V3_18D.resample('60T').mean()

Hg0_SIPEXII = Hg0_SIPEXII.resample('60T').mean()

#------------------------------------------------------------------------------
# RESAMPLE THE BrO DATASETS TO 1-HOUR TIME RESOLUTION

#--------------
# BrO
#--------------
BrO_V1_17   = BrO_V1_17.resample('60T').mean()
BrO_V2_17   = BrO_V2_17.resample('60T').mean()
BrO_V3_17M  = BrO_V3_17M.resample('60T').mean()
BrO_V3_17D  = BrO_V3_17D.resample('60T').mean()

BrO_V1_18   = BrO_V1_18.resample('60T').mean()
BrO_V2_18   = BrO_V2_18.resample('60T').mean()
BrO_V3_18M  = BrO_V3_18M.resample('60T').mean()
BrO_V3_18D  = BrO_V3_18D.resample('60T').mean()

BrO_SIPEXII = BrO_SIPEXII.resample('60T').mean()

#------------------------------------------------------------------------------
# BACK TRAJECTORIES CAMMPCAN (2017-18)

# Set the location for the working directory
os.chdir("/Users/ncp532/Documents/Data/SeaIce_Trajectories/100m/")

# Set the start and end date
DATE_FORMAT = '%Y%m%d%H'
delta_one_hour = timedelta(hours=1)

#-------------
# V1_17
#-------------
# Set the start/end dates
start_V1_17a   = '2017102907'
end_V1_17a     = '2017103123'

start_date    = datetime.strptime(start_V1_17a, DATE_FORMAT)
end_date      = datetime.strptime(end_V1_17a,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1oct0100spring' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V1_17)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V1_17)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V1_17a = pd.concat(Traj)

# Set the start/end dates
start_V1_17b   = '2017110100'
end_V1_17b     = '2017113023'

start_date    = datetime.strptime(start_V1_17b, DATE_FORMAT)
end_date      = datetime.strptime(end_V1_17b,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1nov0100spring' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V1_17)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V1_17)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V1_17b = pd.concat(Traj)

# Set the start/end dates    
start_V1_17c   = '2017120100'
end_V1_17c     = '2017120323'

start_date    = datetime.strptime(start_V1_17c, DATE_FORMAT)
end_date      = datetime.strptime(end_V1_17c,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1dec0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V1_17)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V1_17)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V1_17c = pd.concat(Traj)

# Append V1_17a, V1_17b and V1_17c
Traj_V1_17 = Traj_V1_17a.append(Traj_V1_17b)
Traj_V1_17 = Traj_V1_17.append(Traj_V1_17c)
    
#-------------
# V2_17
#-------------
# Set the start/end dates   
start_V2_17a  = '2017121308'
end_V2_17a    = '2017123123'

start_date    = datetime.strptime(start_V2_17a, DATE_FORMAT)
end_date      = datetime.strptime(end_V2_17a,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1dec0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=8)) # Casey timezone is UT+8
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V2_17)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V2_17)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V2_17a = pd.concat(Traj)

# Set the start/end dates   
start_V2_17b  = '2018010100'
end_V2_17b    = '2018011123'

start_date    = datetime.strptime(start_V2_17b, DATE_FORMAT)
end_date      = datetime.strptime(end_V2_17b,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1jan0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=8)) # Casey timezone is UT+8
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V2_17)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V2_17)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V2_17b = pd.concat(Traj)

# Append V2_17a, V2_17b
Traj_V2_17 = Traj_V2_17a.append(Traj_V2_17b)

#-------------
# V3_17M
#-------------
# Set the start/end dates   
start_V3_17Ma  = '2018011605'
end_V3_17Ma    = '2018013123'

start_date    = datetime.strptime(start_V3_17Ma, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_17Ma,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1jan0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=5)) # Mawson timezone is UT+5
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_17M)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_17M)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)
    
# Combine all the files
Traj_V3_17Ma = pd.concat(Traj)

# Set the start/end dates   
start_V3_17Mb  = '2018020100'
end_V3_17Mb    = '2018022823'

start_date    = datetime.strptime(start_V3_17Mb, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_17Mb,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1feb0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=5)) # Mawson timezone is UT+5
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_17M)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_17M)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)
    
# Combine all the files
Traj_V3_17Mb = pd.concat(Traj)

# Set the start/end dates   
start_V3_17Mc  = '2018030100'
end_V3_17Mc    = '2018030503'

start_date    = datetime.strptime(start_V3_17Mc, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_17Mc,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1mar0100autumn' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=5)) # Mawson timezone is UT+5
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_17M)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_17M)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)
    
# Combine all the files
Traj_V3_17Mc = pd.concat(Traj)

# Append V3_17Ma, V3_17Mb and V3_17Mc
Traj_V3_17M = Traj_V3_17Ma.append(Traj_V3_17Mb)
Traj_V3_17M = Traj_V3_17M.append(Traj_V3_17Mc)

#-------------
# V3_17D
#-------------
# Set the start/end dates   
start_V3_17Da  = '2018011607'
end_V3_17Da    = '2018013123'

start_date    = datetime.strptime(start_V3_17Da, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_17Da,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1jan0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_17D)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_17D)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)
    
# Combine all the files
Traj_V3_17Da = pd.concat(Traj)

# Set the start/end dates   
start_V3_17Db  = '2018020100'
end_V3_17Db    = '2018022823'

start_date    = datetime.strptime(start_V3_17Db, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_17Db,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1feb0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_17D)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_17D)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)
    
# Combine all the files
Traj_V3_17Db = pd.concat(Traj)

# Set the start/end dates   
start_V3_17Dc  = '2018030100'
end_V3_17Dc    = '2018030503'

start_date    = datetime.strptime(start_V3_17Dc, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_17Dc,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1mar0100autumn' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_17D)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_17D)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)
    
# Combine all the files
Traj_V3_17Dc = pd.concat(Traj)

# Append V3_17Da, V3_17Db and V3_17Dc
Traj_V3_17D = Traj_V3_17Da.append(Traj_V3_17Db)
Traj_V3_17D = Traj_V3_17D.append(Traj_V3_17Dc)

#------------------------------------------------------------------------------
# BACK TRAJECTORIES CAMMPCAN (2018-19)

# Set the location for the working directory
os.chdir("/Users/ncp532/Documents/Data/SeaIce_Trajectories/100m/")

#-------------
# V1_18
#-------------
# Set the start/end dates   
start_V1_18a   = '2018102507'
end_V1_18a     = '2018103123'

start_date    = datetime.strptime(start_V1_18a, DATE_FORMAT)
end_date      = datetime.strptime(end_V1_18a,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1oct0100spring' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V1_18)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V1_18)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V1_18a = pd.concat(Traj)

# Set the start/end dates   
start_V1_18b   = '2018110100'
end_V1_18b     = '2018112723'

start_date    = datetime.strptime(start_V1_18b, DATE_FORMAT)
end_date      = datetime.strptime(end_V1_18b,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1nov0100spring' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V1_18)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V1_18)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V1_18b = pd.concat(Traj)

# Append V1_18a and V1_18b
Traj_V1_18 = Traj_V1_18a.append(Traj_V1_18b)

#-------------
# V2_18
#-------------
# Set the start/end dates   
start_V2_18a   = '2018120608'
end_V2_18a     = '2018123123'

start_date    = datetime.strptime(start_V2_18a, DATE_FORMAT)
end_date      = datetime.strptime(end_V2_18a,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1dec0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=8)) # Casey timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V2_18)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V2_18)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V2_18a = pd.concat(Traj)

# Set the start/end dates   
start_V2_18b   = '2019010100'
end_V2_18b     = '2019010723'

start_date    = datetime.strptime(start_V2_18b, DATE_FORMAT)
end_date      = datetime.strptime(end_V2_18b,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1jan0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=8)) # Casey timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V2_18)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V2_18)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V2_18b = pd.concat(Traj)

# Append V2_18a and V2_18b
Traj_V2_18 = Traj_V2_18a.append(Traj_V2_18b)

#-------------
# V3_18M
#-------------
# Set the start/end dates   
start_V3_18Ma  = '2019011305'
end_V3_18Ma    = '2019013123'

start_date    = datetime.strptime(start_V3_18Ma, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_18Ma,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1jan0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=5)) # Mawson timezone is UT+5
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_18M)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_18M)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V3_18Ma = pd.concat(Traj)

# Set the start/end dates   
start_V3_18Mb  = '2019020100'
end_V3_18Mb    = '2019022823'

start_date    = datetime.strptime(start_V3_18Mb, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_18Mb,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1feb0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=5)) # Mawson timezone is UT+5
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_18M)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_18M)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V3_18Mb = pd.concat(Traj)

# Set the start/end dates   
start_V3_18Mc  = '2019030100'
end_V3_18Mc    = '2019030123'

start_date    = datetime.strptime(start_V3_18Mc, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_18Mc,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1mar0100autumn' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=5)) # Mawson timezone is UT+5
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_18M)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_18M)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V3_18Mc = pd.concat(Traj)

# Append V3_18Ma and V3_18Mb
Traj_V3_18M = Traj_V3_18Ma.append(Traj_V3_18Mb)
Traj_V3_18M = Traj_V3_18M.append(Traj_V3_18Mc)

#-------------
# V3_18D
#-------------
# Set the start/end dates   
start_V3_18Da  = '2019011307'
end_V3_18Da    = '2019013123'

start_date    = datetime.strptime(start_V3_18Da, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_18Da,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1jan0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_18D)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_18D)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V3_18Da = pd.concat(Traj)

# Set the start/end dates   
start_V3_18Db  = '2019020100'
end_V3_18Db    = '2019022823'

start_date    = datetime.strptime(start_V3_18Db, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_18Db,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1feb0100summer' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_18D)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_18D)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V3_18Db = pd.concat(Traj)

# Set the start/end dates   
start_V3_18Dc  = '2019030100'
end_V3_18Dc    = '2019030123'

start_date    = datetime.strptime(start_V3_18Dc, DATE_FORMAT)
end_date      = datetime.strptime(end_V3_18Dc,   DATE_FORMAT)

# Save a list of all the file names to the variable all_filenames 
all_filenames = []
date = start_date
while date <= end_date:
    dateA = date.strftime('%Y%m%d%H')
    filename = 'gdas1mar0100autumn' + dateA + '.csv'
    date += delta_one_hour
    all_filenames.append(filename)

# Set an empty array
Traj = []

# Loop over the files in the folder
for f in all_filenames:
    # Set the file name
    file = pd.read_csv(f)
    file['year']     = file['Traj Year']
    file['month']    = file['Traj Mon']
    file['day']      = file['Traj Day']
    file['hour']     = file['Traj Hour']
    file['minute']   = file['Traj Min']
    file['DateTime'] = pd.to_datetime(file[['year', 'month', 'day', 'hour', 'minute']])
    file.index       = (pd.to_datetime(file['DateTime']) + timedelta(hours=7)) # Davis timezone is UT+7
    # Sum the ice contact time (hours)
    file['Over_SeaIce']      = np.sum(file['Traj over Sea Ice?'])
    file['IceContact_100m']  = np.sum(file['Traj over Sea Ice and height < 100 m?'])
    file['IceContact_MLH']   = np.sum(file['Traj over Sea Ice and height < mixed layer top?'])
    # Ice contact time (hours) * sea ice concentration (%)
    file['ContactIcePerc']   = np.sum(file['Traj over Sea Ice?']*file['Sea Ice Conc (0-1)'])
    # Sum the land contact time (hours)
    file['Over_Land']        = np.sum(file['Traj over Ice Sheet?'])
    file['Land_MLH']         = np.sum(file['Traj over Ice Sheet and height < mixed layer top?'])
    # Calculate if traj over ocean
    file['Traj over Ocean?'] = 1 - (file['Traj over Sea Ice?'] + file['Traj over Ice Sheet?'])
    # Sum the ocean contact time (hours)
    file['Over_Ocean']      = np.sum(file['Traj over Ocean?'])
    # Percentage time over Land, Ice, Ocean
    file['Land_Percentage']    = file['Over_Land']/121
    file['SeaIce_Percentage']  = file['Over_SeaIce']/121
    file['Ice100m_Percentage'] = file['IceContact_100m']/121
    file['IceMLH_Percentage']  = file['IceContact_MLH']/121
    file['Ocean_Percentage']   = file['Over_Ocean']/121
    # Weighting for Traj Age
    file['Weighting']        = (file['Traj Age'].values[::-1] - 1) * -1
    # Weighted Land, Ice and Ocean values
    file['Weighted_Ice']     = np.sum(file['Weighting']*file['Traj over Sea Ice?']/121)
    file['Weighted_Ice100m'] = np.sum(file['Weighting']*file['Traj over Sea Ice and height < 100 m?']/121)
    file['Weighted_IceMLH']  = np.sum(file['Weighting']*file['Traj over Sea Ice and height < mixed layer top?']/121)
    file['Weighted_Land']    = np.sum(file['Weighting']*file['Traj over Ice Sheet?']/121)
    file['Weighted_Ocean']   = np.sum(file['Weighting']*file['Traj over Ocean?']/121)
    # Set the Hg0 concentration at the corresponding datetime as a row
    Test = file.join(Hg0_V3_18D)
    file['ng/m3'] = Test.iloc[0]['ng/m3']
    # Set the BrO concentration at the corresponding datetime as a row
    Test = file.join(BrO_V3_18D)
    file['surf_vmr(ppmv)'] = Test.iloc[0]['surf_vmr(ppmv)']*1e6
    # store DataFrame in list
    Traj.append(file)

# Combine all the files
Traj_V3_18Dc = pd.concat(Traj)

# Append V3_18Da and V3_18Db
Traj_V3_18D = Traj_V3_18Da.append(Traj_V3_18Db)
Traj_V3_18D = Traj_V3_18D.append(Traj_V3_18Dc)

#------------------------------------------------------------------------------
# DROP ROWS WHERE BrO CONCENTRATION IS MISSING

# Traj_V1_17  = Traj_V1_17[Traj_V1_17['surf_vmr(ppmv)'].notna()]
# Traj_V2_17  = Traj_V2_17[Traj_V2_17['surf_vmr(ppmv)'].notna()]
# Traj_V3_17M = Traj_V3_17M[Traj_V3_17M['surf_vmr(ppmv)'].notna()]
# Traj_V3_17D = Traj_V3_17D[Traj_V3_17D['surf_vmr(ppmv)'].notna()]

# Traj_V1_18  = Traj_V1_18[Traj_V1_18['surf_vmr(ppmv)'].notna()]
# Traj_V2_18  = Traj_V2_18[Traj_V2_18['surf_vmr(ppmv)'].notna()]
# Traj_V3_18M = Traj_V3_18M[Traj_V3_18M['surf_vmr(ppmv)'].notna()]
# Traj_V3_18D = Traj_V3_18D[Traj_V3_18D['surf_vmr(ppmv)'].notna()]

#------------------------------------------------------------------------------
# GET TRAJECTORIES FOR ONLY TRAJ AGE = 0

Filter1     = Traj_V1_17['Traj Age']  == 0
Traj_V1_17  = Traj_V1_17[Filter1]

Filter2     = Traj_V2_17['Traj Age']  == 0
Traj_V2_17  = Traj_V2_17[Filter2]

Filter3     = Traj_V3_17M['Traj Age'] == 0
Traj_V3_17M = Traj_V3_17M[Filter3]

Filter4     = Traj_V3_17D['Traj Age'] == 0
Traj_V3_17D = Traj_V3_17D[Filter4]

Filter1     = Traj_V1_18['Traj Age']  == 0
Traj_V1_18  = Traj_V1_18[Filter1]

Filter2     = Traj_V2_18['Traj Age']  == 0
Traj_V2_18  = Traj_V2_18[Filter2]

Filter3     = Traj_V3_18M['Traj Age'] == 0
Traj_V3_18M = Traj_V3_18M[Filter3]

Filter4     = Traj_V3_18D['Traj Age'] == 0
Traj_V3_18D = Traj_V3_18D[Filter4]

#------------------------------------------------------------------------------
# REPLACE ERRONEOUS VALUES WITH NAN

Traj_V1_17['Sea Ice Conc (0-1)']  = Traj_V1_17['Sea Ice Conc (0-1)'].replace(999.0,  np.nan)
Traj_V2_17['Sea Ice Conc (0-1)']  = Traj_V2_17['Sea Ice Conc (0-1)'].replace(999.0,  np.nan)
Traj_V3_17M['Sea Ice Conc (0-1)'] = Traj_V3_17M['Sea Ice Conc (0-1)'].replace(999.0, np.nan)
Traj_V3_17D['Sea Ice Conc (0-1)'] = Traj_V3_17D['Sea Ice Conc (0-1)'].replace(999.0, np.nan)

Traj_V1_18['Sea Ice Conc (0-1)']  = Traj_V1_18['Sea Ice Conc (0-1)'].replace(999.0,  np.nan)
Traj_V2_18['Sea Ice Conc (0-1)']  = Traj_V2_18['Sea Ice Conc (0-1)'].replace(999.0,  np.nan)
Traj_V3_18M['Sea Ice Conc (0-1)'] = Traj_V3_18M['Sea Ice Conc (0-1)'].replace(999.0, np.nan)
Traj_V3_18D['Sea Ice Conc (0-1)'] = Traj_V3_18D['Sea Ice Conc (0-1)'].replace(999.0, np.nan)

#------------------------------------------------------------------------------
# EXPORT AS .CSV FILE

# Save to .csv
Traj_V1_17.to_csv('/Users/ncp532/Desktop/Traj_V1_17.csv')
Traj_V2_17.to_csv('/Users/ncp532/Desktop/Traj_V2_17.csv')
Traj_V3_17M.to_csv('/Users/ncp532/Desktop/Traj_V3_17M.csv')
Traj_V3_17D.to_csv('/Users/ncp532/Desktop/Traj_V3_17D.csv')

Traj_V1_18.to_csv('/Users/ncp532/Desktop/Traj_V1_18.csv')
Traj_V2_18.to_csv('/Users/ncp532/Desktop/Traj_V2_18.csv')
Traj_V3_18M.to_csv('/Users/ncp532/Desktop/Traj_V3_18M.csv')
Traj_V3_18D.to_csv('/Users/ncp532/Desktop/Traj_V3_18D.csv')

#------------------------------------------------------------------------------
# PLOT THE GRAPH
fig = plt.figure()

gs = gridspec.GridSpec(nrows=2,
                       ncols=1, 
                       figure=fig)

#-----------------------------
# Graph 1
ax1 = plt.subplot(gs[0,0])
ax2 = ax1.twinx()
ax2.spines["right"].set_color('blue')
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.04))
ax3.spines["right"].set_color('magenta')
ax4 = ax1.twinx()
ax4.spines["right"].set_position(("axes", 1.08))
ax4.spines["right"].set_color('grey')

# Plot the variables
ax1.plot(Traj_V3_18M.index, Traj_V3_18M['ng/m3'], marker='None', ls='-', c='black', label ="V3 (2018-19)")
ax2.plot(Traj_V3_18M.index, Traj_V3_18M['surf_vmr(ppmv)'], marker='None', ls='-', c='blue', label ="V3 (2018-19)")
ax3.plot(Traj_V3_18M.index, Traj_V3_18M['Traj Lat'], ls='--', c='magenta', label ='Latitude')
ax4.plot(Traj_V3_18M.index, Traj_V3_18M['Sea Ice Conc (0-1)']*100, ls='--', c='grey', label ='Sea Ice Concentration')

