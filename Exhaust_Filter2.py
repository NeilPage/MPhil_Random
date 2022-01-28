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

# Exhaust Filter
exhaust = xr.open_dataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/Exhaust_ID/AAS_4292_ExhaustID_201718_AA_MARCUS.nc')
exhaust = exhaust.to_dataframe()

# O3 Data
V1_O3 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/All_O3_1sec.csv', index_col=0)

# CO Data
V1_CO = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/All_CO_1sec.csv', index_col=0)

#------------------------------------------------------------------------------
# SET THE DATE

exhaust.index   = pd.to_datetime(exhaust.index,     dayfirst=True)
V1_O3.index     = pd.to_datetime(V1_O3.index,       dayfirst=True)
V1_CO.index     = pd.to_datetime(V1_CO.index,       dayfirst=True)

#------------------------------------------------------------------------------
# CLEAN UP THE O3 DATA (REMOVE ERRONEOUS DATA)

# O3 (positive values only)
filter1 = V1_O3['O3_(ppb)'] >= 0
V1_O3   = V1_O3[filter1]

# O3 (get rid of stupidly high values)
filter2 = V1_O3['O3_(ppb)'] <= 50
V1_O3   = V1_O3[filter2]

# O3 (drop nan values)
V1_O3   = V1_O3[V1_O3['O3_(ppb)'].notna()]

#------------------------------------------------------------------------------
# MERGE THE DATASETS WITH THE EXHAUST FILTER

#Filtered_NO2 = V1_NO2.merge(exhaust, how='left',left_index=True, right_index=True)
Filtered_O3 = V1_O3.merge(exhaust, how='left',left_index=True, right_index=True)
Filtered_CO = V1_CO.merge(exhaust, how='left',left_index=True, right_index=True)

#------------------------------------------------------------------------------
# FILTER FOR THE EXHAUST (exhaust_4mad02thresh)

# O3
is_false = Filtered_O3['exhaust_4mad02thresh']==False
O3a      = Filtered_O3[is_false]

# CO
is_false = Filtered_CO['exhaust_4mad02thresh']==False
COa      = Filtered_CO[is_false]

#------------------------------------------------------------------------------
# REAVERAGE THE DATASETS FROM 1 SECOND TO 60 MIN

O3b    = O3a.resample('1T').mean()
COb    = COa.resample('1T').mean()

V1_O3b = V1_O3.resample('1T').mean()
V1_COb = V1_CO.resample('1T').mean()

#------------------------------------------------------------------------------
# SEPERATE THE DATASETS INTO THE SEPERATE VOYAGES

#-----------------------------
# V1_17 Davis (29 Oct - 3 Dec 2017)
#-----------------------------
start_date = '2017-10-29'
end_date   = '2017-12-04'

Davis      = (O3b.index >= start_date) & (O3b.index < end_date)
O3_V1_17   = O3b[Davis]

#-----------------------------
# V2_17 Casey (13 Dec - 11 Jan 2018)
#-----------------------------
start_date = '2017-12-13'
end_date   = '2018-01-12'

Casey      = (O3b.index >= start_date) & (O3b.index < end_date)
O3_V2_17   = O3b[Casey]

#-----------------------------
# V3_17 Mawson & Davis (16 Jan - 7 Mar 2018)
#-----------------------------
start_date = '2018-01-16'
end_date   = '2018-03-07'

Mawson     = (O3b.index >= start_date) & (O3b.index < end_date)
O3_V3_17   = O3b[Mawson]

#------------------------------------------------------------------------------
# SAVE THE FILTERED DATAFRAME AS .CSV

O3_V1_17.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/V1_O3_1min_Fil.csv')
O3_V2_17.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/V2_O3_1min_Fil.csv')
O3_V3_17.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/V3_O3_1min_Fil.csv')

#------------------------------------------------------------------------------
# PLOT THE GRAPH
fig = plt.figure()

gs = gridspec.GridSpec(nrows=2,
                       ncols=1, 
                       figure=fig)

#-----------------------------
# Graph 1
ax1 = plt.subplot(gs[0,0])

# Plot the variables
ax1.plot(V1_O3b.index, V1_O3b['O3_(ppb)'], marker='o', ls='None', ms=2, c='black', label ="O3 (All")
ax1.plot(O3b.index,    O3b['O3_(ppb)'],    marker='o', ls='None', ms=2, c='red',   label ="O3 (Filtered")

#-----------------------------
# Graph 2
ax1 = plt.subplot(gs[1,0])

# Plot the variables
ax1.plot(V1_COb.index, V1_COb['CO'], marker='o', ls='None', ms=2, c='black', label ="CO (All")
ax1.plot(COb.index,    COb['CO'],    marker='o', ls='None', ms=2, c='red',   label ="CO (Filtered")
ax1.set_yscale('log')

