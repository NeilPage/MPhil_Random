#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:03:02 2020

@author: ncp532
"""

# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from matplotlib import gridspec

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
import datetime as dt
from datetime import datetime,time, timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# Aerosol Size Distribution

#---------
# BrO
#---------
# BrO VMR
V1_2017 = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V1_17/all_BrO/V1_17_BrO_VMR.csv', index_col=0) # BrO data for CAMPCANN V1 (2017/18)

# BrO Error
V1_Error = pd.read_csv('/Users/ncp532/Documents/data/V1_17_APriori/V1_17/BrO_error/V1_17_BrO_error.csv', index_col=0) # BrO error data for CAMPCANN V1 (2017/18)

# Calculate the Relative Error (>=0.6)
Filter = V1_Error / V1_2017

# Apply the filter
V1_17F = Filter < 0.6
V1_17T = V1_2017[V1_17F]

#---------
# Aerosol Number Concentration
#---------
DF1 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/AEROSOL_SIZE_DISTRIBUTION/ASD/Vovage_1_WIBS_nc.csv')

#---------
# Sea Salt Aerosol
#---------
DF2 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/AEROSOL_SIZE_DISTRIBUTION/ASD/acsm_clean.csv')
        
#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

#------------------
# BrO
y = V1_17T.index # set the values for the y-axis
x = np.array(V1_17T.dtypes.index) # set the values for the x-axis
z = V1_17T.copy() # identify the matrix containing the z-values (BrO in ppMv)
z[z==-9999]=np.nan # set the erroneous values as NaN 
z = z.loc[:]*1e6 # change from ppMv to pptv
mz=np.ma.masked_where(np.isnan(z),z) 

# Date
Date1 = V1_2017.columns.values
date1=[]
for i in range(len(Date1)):
    date1.append(datetime.strptime(Date1[i],'%d/%m/%Y %H:%M:%S')+ timedelta(hours=7))

#------------------
# Aerosol Number Concentration

# Date
Date2 = np.array(DF1['DateTime'])
date2=[]
for i in range(len(Date2)):
    date2.append(datetime.strptime(Date2[i],'%d/%m/%Y %H:%M:%S')+ timedelta(hours=7))

# Aero
Aero = np.array(DF1['All'])

#------------------
# Sea Salt Aerosol

# Date
Date3 = np.array(DF2['DateTime'])
date3=[]
for i in range(len(Date3)):
    date3.append(datetime.strptime(Date3[i],'%d/%m/%Y %H:%M:%S')+ timedelta(hours=7))

# Aero
SSA = np.array(DF2['SSA'])

#------------------------------------------------------------------------------
# PLOT THE GRAPH

fig = plt.figure(figsize=(10,6))
ax  = fig.add_subplot(111)
ax2 = ax.twinx()

cmap=plt.cm.jet
norm = BoundaryNorm(np.arange(0,11,1), cmap.N)

col1 = ax.pcolormesh(date1, y, mz, vmin=0, vmax=10, norm=norm, cmap=cmap)
ax2.scatter(date2, Aero, marker='x', color='magenta', label ='Aerosol')

plt.xlim(datetime(2017,11,14,0,0,0),datetime(2017,11,23,23,59,59))

ax2.set_ylabel('Altitude (km)', fontsize=20)
ax2.set_ylabel('Number Concentration\n(0.5 to 10 \u03BCm)', fontsize=20)
ax.set_xlabel('Date', fontsize=20)

# Format y-axis
ax.set_ylim(0.1,3.0) # On Station
ax2.set_ylim(0,50)
fig.autofmt_xdate()

# Format x-axis
xmajor_formatter = mdates.DateFormatter('%d %b') # format how the date is displayed
ax.xaxis.set_major_formatter(xmajor_formatter)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=2)) # set the interval between dispalyed dates
ax.xaxis.set_minor_locator(mdates.HourLocator(interval=4))

