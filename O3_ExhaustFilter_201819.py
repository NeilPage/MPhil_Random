#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:59:25 2019

@author: ncp532
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:16:43 2019

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files
import xarray as xr                 # xarray to open multiple netcdf files    

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
# DEFINE THE DATASET

# O3
O3_V1_18   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/O3/V1_O3_1min.csv', index_col=0)
O3_V1_18.rename(columns={'O3':'O3_(ppb)'},inplace=True) # rename the column from 'O3' to 'O3_(ppb)'
O3_V2_18   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/O3/V2_O3_1min.csv', index_col=0)
O3_V2_18.rename(columns={'O3':'O3_(ppb)'},inplace=True) # rename the column from 'O3' to 'O3_(ppb)'
O3_V2_18   = O3_V2_18.loc[~O3_V2_18.index.duplicated(keep='first')] # remove duplicate values from the .csv file
O3_V3_18   = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/O3/V3_O3_1min.csv', index_col=0)
O3_V3_18   = O3_V3_18.loc[~O3_V3_18.index.duplicated(keep='first')] # remove duplicate values from the .csv file
O3_V3_18.rename(columns={'O3':'O3_(ppb)'},inplace=True) # rename the column from 'O3' to 'O3_(ppb)'

# Met
Met_V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/ShipTrack/V1_18_underway_60.csv', index_col=0) 
Met_V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/ShipTrack/V2_18_underway_60.csv', index_col=0)
Met_V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/ShipTrack/V3_18_underway_60.csv', index_col=0) 

#------------------------------------------------------------------------------
# SET THE DATE

# O3
O3_V1_18.index  = pd.to_datetime(O3_V1_18.index,  dayfirst=True)
O3_V2_18.index  = pd.to_datetime(O3_V2_18.index,  dayfirst=True)
O3_V3_18.index  = pd.to_datetime(O3_V3_18.index,  dayfirst=True)

# Met
Met_V1_18.index = pd.to_datetime(Met_V1_18.index, dayfirst=True)
Met_V2_18.index = pd.to_datetime(Met_V2_18.index, dayfirst=True)
Met_V3_18.index = pd.to_datetime(Met_V3_18.index, dayfirst=True)

#------------------------------------------------------------------------------
# RESAMPLE THE MET DATASETS TO 1 MIN

Met_V1_18 = Met_V1_18.resample('1T').mean()
Met_V2_18 = Met_V2_18.resample('1T').mean()
Met_V3_18 = Met_V3_18.resample('1T').mean()

#------------------------------------------------------------------------------
# Add Wind Speed Average & Wind Direction Average

# Wind Speed
Met_V1_18['WS_Avg'] = Met_V1_18[['wnd_spd_port_corr_knot', 'wnd_spd_strbrd_corr_knot']].mean(axis=1) * 0.514444444 # Met wind speed average (port & strbrd side) (m/s)
Met_V2_18['WS_Avg'] = Met_V2_18[['wnd_spd_port_corr_knot', 'wnd_spd_strbrd_corr_knot']].mean(axis=1) * 0.514444444 # Met wind speed average (port & strbrd side) (m/s)
Met_V3_18['WS_Avg'] = Met_V3_18[['wnd_spd_port_corr_knot', 'wnd_spd_strbrd_corr_knot']].mean(axis=1) * 0.514444444 # Met wind speed average (port & strbrd side) (m/s)

# Wind Direction
Met_V1_18['WD_Avg'] = Met_V1_18[['wnd_dir_port_corr_deg',  'wnd_dir_strbrd_corr_deg']].mean(axis=1) # Wind direction (degrees)
Met_V2_18['WD_Avg'] = Met_V2_18[['wnd_dir_port_corr_deg',  'wnd_dir_strbrd_corr_deg']].mean(axis=1) # Wind direction (degrees)
Met_V3_18['WD_Avg'] = Met_V3_18[['wnd_dir_port_corr_deg',  'wnd_dir_strbrd_corr_deg']].mean(axis=1) # Wind direction (degrees)

#------------------------------------------------------------------------------
# MERGE THE O3 & Met DATASETS

O3_V1_18 = pd.concat([O3_V1_18, Met_V1_18],   axis=1, join='inner')
O3_V2_18 = pd.concat([O3_V2_18, Met_V2_18],   axis=1, join='inner')
O3_V3_18 = pd.concat([O3_V3_18, Met_V3_18],   axis=1, join='inner')

#------------------------------------------------------------------------------
# FILTER FOR THE SHIPS EXHAUST

# O3 (remove data when wind direction is 90-270 degrees)
O3_90     = O3_V1_18['WD_Avg'] <=90 # <=90
D10       = O3_V1_18[O3_90]
O3_270    = O3_V1_18['WD_Avg'] >=270 # >=270 
D11       = O3_V1_18[O3_270]
O3_V1_18b = pd.concat([D10,D11], axis=0)

O3_90     = O3_V2_18['WD_Avg'] <=90 # <=90
D10       = O3_V2_18[O3_90]
O3_270    = O3_V2_18['WD_Avg'] >=270 # >=270 
D11       = O3_V2_18[O3_270]
O3_V2_18b = pd.concat([D10,D11], axis=0)

O3_90     = O3_V3_18['WD_Avg'] <=90 # <=90
D10       = O3_V3_18[O3_90]
O3_270    = O3_V3_18['WD_Avg'] >=270 # >=270 
D11       = O3_V3_18[O3_270]
O3_V3_18b = pd.concat([D10,D11], axis=0)

# Step 3 (Remove data when wind speed is below 5 knots (or 2.57222222 m/s))
O3_5knot  = O3_V1_18b['WS_Avg'] >=2.57222222
O3_V1_18c = O3_V1_18b[O3_5knot]

O3_5knot  = O3_V2_18b['WS_Avg'] >=2.57222222
O3_V2_18c = O3_V2_18b[O3_5knot]

O3_5knot  = O3_V3_18b['WS_Avg'] >=2.57222222
O3_V3_18c = O3_V3_18b[O3_5knot]

#------------------------------------------------------------------------------
# REAVERAGE THE DATASETS FROM 1 MIN TO 60 MIN

V1_O3 = O3_V1_18b.resample('60T').mean()
V2_O3 = O3_V2_18b.resample('60T').mean()
V3_O3 = O3_V3_18b.resample('60T').mean()

#------------------------------------------------------------------------------
# SAVE THE FILTERED DATAFRAME AS .CSV

O3_V1_18b.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/O3/V1_O3_1min_Fil.csv')
O3_V2_18b.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/O3/V2_O3_1min_Fil.csv')
O3_V3_18b.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/O3/V3_O3_1min_Fil.csv')

#------------------------------------------------------------------------------
# PLOT THE GRAPH
fig = plt.figure()

gs = gridspec.GridSpec(nrows=3,
                       ncols=1, 
                       figure=fig)

#-----------------------------
# Graph 1
ax1 = plt.subplot(gs[0,0])

# Plot the variables
ax1.plot(O3_V1_18.index,  O3_V1_18['O3_(ppb)'],  marker='o', ls='None', ms=2, c='black', label ="O3 (All")
ax1.plot(V1_O3.index,     V1_O3['O3_(ppb)'],     marker='o', ls='None', ms=2, c='red',   label ="O3 (All")
#ax1.plot(O3_V1_18b.index, O3_V1_18b['O3_(ppb)'], marker='o', ls='None', ms=2, c='red',   label ="O3 (Filtered")
#ax1.plot(O3_V1_18c.index, O3_V1_18c['O3_(ppb)'], marker='o', ls='None', ms=2, c='blue',  label ="O3 (Filtered")

#-----------------------------
# Graph 2
ax1 = plt.subplot(gs[1,0])

# Plot the variables
ax1.plot(O3_V2_18.index,  O3_V2_18['O3_(ppb)'],  marker='o', ls='None', ms=2, c='black', label ="O3 (All")
ax1.plot(V2_O3.index,     V2_O3['O3_(ppb)'],     marker='o', ls='None', ms=2, c='red',   label ="O3 (All")
#ax1.plot(O3_V2_18b.index, O3_V2_18b['O3_(ppb)'], marker='o', ls='None', ms=2, c='red',   label ="O3 (Filtered")
#ax1.plot(O3_V2_18c.index, O3_V2_18c['O3_(ppb)'], marker='o', ls='None', ms=2, c='blue',  label ="O3 (Filtered")

#-----------------------------
# Graph 3
ax1 = plt.subplot(gs[2,0])

# Plot the variables
ax1.plot(O3_V3_18.index,  O3_V3_18['O3_(ppb)'],  marker='o', ls='None', ms=2, c='black', label ="O3 (All")
ax1.plot(V3_O3.index,     V3_O3['O3_(ppb)'],     marker='o', ls='None', ms=2, c='red',   label ="O3 (All")
#ax1.plot(O3_V3_18b.index, O3_V3_18b['O3_(ppb)'], marker='o', ls='None', ms=2, c='red',   label ="O3 (Filtered")
#ax1.plot(O3_V3_18c.index, O3_V3_18c['O3_(ppb)'], marker='o', ls='None', ms=2, c='blue',  label ="O3 (Filtered")

