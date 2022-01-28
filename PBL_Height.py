#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:20:49 2020

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files
import xarray as xr

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# DEFINE THE DATASETS

#--------------
# MERRA-2
#--------------

MERRA2_V1_17   = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V1_17_MERRA2.csv',   index_col=0)
MERRA2_V2_17   = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V2_17_MERRA2.csv',   index_col=0)
MERRA2_V3_17M  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V3_17M_MERRA2.csv',  index_col=0)
MERRA2_V3_17D  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V3_17D_MERRA2.csv',  index_col=0)

MERRA2_V1_18   = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V1_18_MERRA2.csv',   index_col=0) 
MERRA2_V2_18   = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V2_18_MERRA2.csv',   index_col=0)
MERRA2_V3_18M  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V3_18M_MERRA2.csv',  index_col=0) 
MERRA2_V3_18D  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/V3_18D_MERRA2.csv',  index_col=0)

MERRA2_SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/SIPEXII_MERRA2.csv', index_col=0)

#--------------
# Radiosonde
#--------------
RS_V1_17  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde.csv', index_col=0)
RS_V2_17  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde.csv', index_col=0)
RS_V3_17M = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde.csv', index_col=0)
RS_V3_17D = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde.csv', index_col=0)

RS_Interp_V1_17  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde_Interp1hr.csv', index_col=0)
RS_Interp_V2_17  = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde_Interp1hr.csv', index_col=0)
RS_Interp_V3_17M = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde_Interp1hr.csv', index_col=0)
RS_Interp_V3_17D = pd.read_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde_Interp1hr.csv', index_col=0)

#------------------------------------------------------------------------------
# SET THE DATE

#--------------
# MERRA-2
#--------------

MERRA2_V1_17.index   = (pd.to_datetime(MERRA2_V1_17.index,   dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7
MERRA2_V2_17.index   = (pd.to_datetime(MERRA2_V2_17.index,   dayfirst=True) + timedelta(hours=8)) # Casey timezone is UT+8
MERRA2_V3_17M.index  = (pd.to_datetime(MERRA2_V3_17M.index,  dayfirst=True) + timedelta(hours=5)) # Mawson timezone is UT+5
MERRA2_V3_17D.index  = (pd.to_datetime(MERRA2_V3_17D.index,  dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7

MERRA2_V1_18.index   = (pd.to_datetime(MERRA2_V1_18.index,   dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7
MERRA2_V2_18.index   = (pd.to_datetime(MERRA2_V2_18.index,   dayfirst=True) + timedelta(hours=8)) # Casey timezone is UT+8
MERRA2_V3_18M.index  = (pd.to_datetime(MERRA2_V3_18M.index,  dayfirst=True) + timedelta(hours=5)) # Mawson timezone is UT+5
MERRA2_V3_18D.index  = (pd.to_datetime(MERRA2_V3_18D.index,  dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7

MERRA2_SIPEXII.index = (pd.to_datetime(MERRA2_SIPEXII.index, dayfirst=True) + timedelta(hours=8)) # SIPEXII timezone is UT+8

#--------------
# Radiosonde
#--------------
RS_V1_17.index   = (pd.to_datetime(RS_V1_17.index,   dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7
RS_V2_17.index   = (pd.to_datetime(RS_V2_17.index,   dayfirst=True) + timedelta(hours=8)) # Casey timezone is UT+8
RS_V3_17M.index  = (pd.to_datetime(RS_V3_17M.index,  dayfirst=True) + timedelta(hours=5)) # Mawson timezone is UT+5
RS_V3_17D.index  = (pd.to_datetime(RS_V3_17D.index,  dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7

RS_Interp_V1_17.index   = (pd.to_datetime(RS_Interp_V1_17.index,   dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7
RS_Interp_V2_17.index   = (pd.to_datetime(RS_Interp_V2_17.index,   dayfirst=True) + timedelta(hours=8)) # Casey timezone is UT+8
RS_Interp_V3_17M.index  = (pd.to_datetime(RS_Interp_V3_17M.index,  dayfirst=True) + timedelta(hours=5)) # Mawson timezone is UT+5
RS_Interp_V3_17D.index  = (pd.to_datetime(RS_Interp_V3_17D.index,  dayfirst=True) + timedelta(hours=7)) # Davis timezone is UT+7

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

#--------------
# MERRA-2
#--------------
MLH_V1_17   = MERRA2_V1_17['MLH']*1000   # Richardson mixed layer height (m)
MLH_V2_17   = MERRA2_V2_17['MLH']*1000   # Richardson mixed layer height (m)
MLH_V3_17M  = MERRA2_V3_17M['MLH']*1000  # Richardson mixed layer height (m)
MLH_V3_17D  = MERRA2_V3_17D['MLH']*1000  # Richardson mixed layer height (m)

MLH_V1_18   = MERRA2_V1_18['MLH']*1000   # Richardson mixed layer height (m)
MLH_V2_18   = MERRA2_V2_18['MLH']*1000   # Richardson mixed layer height (m)
MLH_V3_18M  = MERRA2_V3_18M['MLH']*1000  # Richardson mixed layer height (m)
MLH_V3_18D  = MERRA2_V3_18D['MLH']*1000  # Richardson mixed layer height (m)

MLH_SIPEXII = MERRA2_SIPEXII['MLH']*1000 # Richardson mixed layer height (m)

#--------------
# Radiosonde
#--------------
MLH_V1_17_RS  = RS_V1_17['MLH']*1000   # Richardson mixed layer height (m)
MLH_V2_17_RS  = RS_V2_17['MLH']*1000   # Richardson mixed layer height (m)
MLH_V3_17M_RS = RS_V3_17M['MLH']*1000  # Richardson mixed layer height (m)
MLH_V3_17D_RS = RS_V3_17D['MLH']*1000  # Richardson mixed layer height (m)

MLH_V1_17_RS_Interp  = RS_Interp_V1_17['MLH']*1000   # Richardson mixed layer height (m)
MLH_V2_17_RS_Interp  = RS_Interp_V2_17['MLH']*1000   # Richardson mixed layer height (m)
MLH_V3_17M_RS_Interp = RS_Interp_V3_17M['MLH']*1000  # Richardson mixed layer height (m)
MLH_V3_17D_RS_Interp = RS_Interp_V3_17D['MLH']*1000  # Richardson mixed layer height (m)

#------------------------------------------------------------------------------
# PLOT THE GRAPH

fig = plt.figure()

# gs = gridspec.GridSpec(nrows=4,
#                        ncols=3, 
#                        figure=fig, 
#                        width_ratios= [0.25,0.25,0.25,0.25],
#                        height_ratios=[0.33, 0.33, 0.33])

#-----------------------------
# Graph 1
#ax1 = plt.subplot(gs[0,0:3])
ax1 = plt.subplot(311)

# Plot the variables
ax1.plot(MERRA2_V1_17.index, MLH_V1_17, marker='x', markersize = 2.0, ls='None', c='blue', label ='MERRA-2')
ax1.scatter(RS_V1_17.index, MLH_V1_17_RS, marker='x', ls='None', c='orange', label ='RadioSonde')

# Format x-axis
plt.xlim(datetime(2017,11,1,0,0,0),datetime(2017,12,13,23,59,59)) # all dates
xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
ax1.xaxis.set_major_formatter(xmajor_formatter)
xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
ax1.xaxis.set_minor_formatter(xminor_formatter)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=7)) # set the interval between dispalyed dates
ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax1.tick_params(axis='x',pad=15)

# Format y-axis
ax1.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
ax1.yaxis.label.set_color('black')
ax1.tick_params(axis='y', which='both', colors='black')
ax1.set_ylim(0,)

# Plot the axis labels
ax1.set_ylabel('Altitude (m)', fontsize=15)
#ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
#plt.legend(bbox_to_anchor=(1.01, 0.95), loc=2, borderaxespad=0., title='V1 (2017/18)')
plt.title('Richardson mixing layer height', fontsize=25, y=1.05)#, x = 1.05)

# Text box in upper left
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.025, 0.925, "V1 (2017-18)", transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)

#-----------------------------
# Graph 2
#ax1 = plt.subplot(gs[1,0:4])
ax1 = plt.subplot(312)

# Plot the variables
ax1.plot(MERRA2_V2_17.index, MLH_V2_17, marker='x', markersize = 2.0, ls='None', c='blue', label ='MERRA-2')
ax1.scatter(RS_V2_17.index, MLH_V2_17_RS, marker='x', ls='None', c='orange', label ='RadioSonde')

# Format x-axis
plt.xlim(datetime(2017,12,13,0,0,0),datetime(2018,1,16,23,59,59)) # all dates
xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
ax1.xaxis.set_major_formatter(xmajor_formatter)
xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
ax1.xaxis.set_minor_formatter(xminor_formatter)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=7)) # set the interval between dispalyed dates
ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax1.tick_params(axis='x',pad=15)

# Format y-axis
ax1.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
ax1.yaxis.label.set_color('black')
ax1.tick_params(axis='y', which='both', colors='black')
ax1.set_ylim(0,)

# Plot the axis labels
ax1.set_ylabel('Altitude (m)', fontsize=15)
#ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)

# Text box in upper left
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.025, 0.925, "V2 (2017-18)", transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)

#-----------------------------
# Graph 3
#ax1 = plt.subplot(gs[2,0:4])
ax1 = plt.subplot(313)

# Plot the variables
ax1.plot(MERRA2_V3_17M.index, MLH_V3_17M, marker='x', markersize = 2.0, ls='None', c='blue', label ='MERRA-2')
ax1.scatter(RS_V3_17M.index, MLH_V3_17M_RS, marker='x', ls='None', c='orange', label ='RadioSonde')

# Format x-axis
plt.xlim(datetime(2018,1,16,0,0,0),datetime(2018,3,1,23,59,59)) # all dates
xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
ax1.xaxis.set_major_formatter(xmajor_formatter)
xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
ax1.xaxis.set_minor_formatter(xminor_formatter)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=7)) # set the interval between dispalyed dates
ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax1.tick_params(axis='x',pad=15)

# Format y-axis
ax1.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
ax1.yaxis.label.set_color('black')
ax1.tick_params(axis='y', which='both', colors='black')
ax1.set_ylim(0,)

# Plot the axis labels
ax1.set_ylabel('Altitude (m)', fontsize=15)
ax1.set_xlabel('Date', fontsize=15, labelpad=-0.1)
plt.legend(bbox_to_anchor=(1.01, 0.95), loc=2, borderaxespad=0., fontsize=15)

# Text box in upper left
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.025, 0.925, "V3 (2017-18)", transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)


# #-----------------------------
# # Graph 4
# ax1 = plt.subplot(gs[3,0:4])

# # Plot the variables
# ax1.plot(MERRA2_V3_17D.index, MLH_V3_17D, marker='x', markersize = 1.0, ls='None', c='blue', label ='MERRA-2')
# ax1.scatter(RS_V3_17M.index, MLH_V3_17M_RS, marker='x', ls='None', c='orange', label ='RadioSonde')

# # Format x-axis
# plt.xlim(datetime(2018,1,16,0,0,0),datetime(2018,3,1,23,59,59)) # all dates
# xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
# ax1.xaxis.set_major_formatter(xmajor_formatter)
# xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
# ax1.xaxis.set_minor_formatter(xminor_formatter)
# ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
# ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
# ax1.tick_params(axis='x',pad=15)

# # Format y-axis
# ax1.yaxis.set_major_locator(ticker.MultipleLocator(500))
# ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
# ax1.yaxis.label.set_color('black')
# ax1.tick_params(axis='y', which='both', colors='black')
# ax1.set_ylim(0,)

# # Plot the axis labels
# ax1.set_ylabel('ALtitude (m)', fontsize=12)
# #ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
# plt.legend(bbox_to_anchor=(1.01, 0.95), loc=2, borderaxespad=0., title='V3D (2017/18)')

# # #-----------------------------
# # # Graph 5
# # ax1 = plt.subplot(gs[0,2:4])

# # # Plot the variables
# # ax1.plot(MERRA2_V1_18.index, MLH_V1_18, marker='x', markersize = 1.0, ls='None', c='red', label ='V1 (2018/19)')

# # # Format x-axis
# # #plt.xlim(datetime(2018,1,2),datetime(2018,1,3)) # all dates
# # xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
# # ax1.xaxis.set_major_formatter(xmajor_formatter)
# # xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
# # ax1.xaxis.set_minor_formatter(xminor_formatter)
# # ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
# # ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
# # ax1.tick_params(axis='x',pad=15)

# # # Format y-axis
# # ax1.yaxis.set_major_locator(ticker.MultipleLocator(200))
# # ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
# # ax1.yaxis.label.set_color('black')
# # ax1.tick_params(axis='y', which='both', colors='black')
# # ax1.set_ylim(0,)

# # # Plot the axis labels
# # #ax1.set_ylabel('ALtitude (m)', fontsize=12)
# # #ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
# # plt.legend(bbox_to_anchor=(0.70, 0.95), loc=2, borderaxespad=0.)

# # #-----------------------------
# # # Graph 6
# # ax1 = plt.subplot(gs[1,2:4])

# # # Plot the variables
# # ax1.plot(MERRA2_V2_18.index, MLH_V2_18, marker='x', markersize = 1.0, ls='None', c='red', label ='V2 (2018/19)')

# # # Format x-axis
# # #plt.xlim(datetime(2018,1,2),datetime(2018,1,3)) # all dates
# # xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
# # ax1.xaxis.set_major_formatter(xmajor_formatter)
# # xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
# # ax1.xaxis.set_minor_formatter(xminor_formatter)
# # ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
# # ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
# # ax1.tick_params(axis='x',pad=15)

# # # Format y-axis
# # ax1.yaxis.set_major_locator(ticker.MultipleLocator(100))
# # ax1.yaxis.set_minor_locator(ticker.MultipleLocator(25))
# # ax1.yaxis.label.set_color('black')
# # ax1.tick_params(axis='y', which='both', colors='black')
# # ax1.set_ylim(0,)

# # # Plot the axis labels
# # #ax1.set_ylabel('ALtitude (m)', fontsize=12)
# # #ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
# # plt.legend(bbox_to_anchor=(0.25, 0.95), loc=2, borderaxespad=0.)

# # #-----------------------------
# # # Graph 7
# # ax1 = plt.subplot(gs[2,2:4])

# # # Plot the variables
# # ax1.plot(MERRA2_V3_18M.index, MLH_V3_18M, marker='x', markersize = 1.0, ls='None', c='red', label ='V3M (2018/19)')

# # # Format x-axis
# # #plt.xlim(datetime(2018,1,2),datetime(2018,1,3)) # all dates
# # xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
# # ax1.xaxis.set_major_formatter(xmajor_formatter)
# # xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
# # ax1.xaxis.set_minor_formatter(xminor_formatter)
# # ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
# # ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
# # ax1.tick_params(axis='x',pad=15)

# # # Format y-axis
# # ax1.yaxis.set_major_locator(ticker.MultipleLocator(200))
# # ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
# # ax1.yaxis.label.set_color('black')
# # ax1.tick_params(axis='y', which='both', colors='black')
# # ax1.set_ylim(0,)

# # # Plot the axis labels
# # #ax1.set_ylabel('ALtitude (m)', fontsize=12)
# # #ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
# # plt.legend(bbox_to_anchor=(0.25, 0.95), loc=2, borderaxespad=0.)

# # #-----------------------------
# # # Graph 8
# # ax1 = plt.subplot(gs[3,2:4])

# # # Plot the variables
# # ax1.plot(MERRA2_V3_18D.index, MLH_V3_18D, marker='x', markersize = 1.0, ls='None', c='red', label ='V3D (2018/19)')

# # # Format x-axis
# # #plt.xlim(datetime(2018,1,2),datetime(2018,1,3)) # all dates
# # xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
# # ax1.xaxis.set_major_formatter(xmajor_formatter)
# # xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
# # ax1.xaxis.set_minor_formatter(xminor_formatter)
# # ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
# # ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
# # ax1.tick_params(axis='x',pad=15)

# # # Format y-axis
# # ax1.yaxis.set_major_locator(ticker.MultipleLocator(100))
# # ax1.yaxis.set_minor_locator(ticker.MultipleLocator(25))
# # ax1.yaxis.label.set_color('black')
# # ax1.tick_params(axis='y', which='both', colors='black')
# # ax1.set_ylim(0,)

# # # Plot the axis labels
# # #ax1.set_ylabel('ALtitude (m)', fontsize=12)
# # #ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
# # plt.legend(bbox_to_anchor=(0.40, 0.95), loc=2, borderaxespad=0.)

# # #-----------------------------
# # # Graph 9
# # ax1 = plt.subplot(gs[4,1:3])

# # # Plot the variables
# # ax1.plot(MERRA2_SIPEXII.index, MLH_SIPEXII, marker='x', markersize = 1.0, ls='None', c='green', label ='SIPEXII (2012)')

# # # Format x-axis
# # #plt.xlim(datetime(2018,1,2),datetime(2018,1,3)) # all dates
# # xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
# # ax1.xaxis.set_major_formatter(xmajor_formatter)
# # xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
# # ax1.xaxis.set_minor_formatter(xminor_formatter)
# # ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
# # ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
# # ax1.tick_params(axis='x',pad=15)

# # # Format y-axis
# # ax1.yaxis.set_major_locator(ticker.MultipleLocator(300))
# # ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
# # ax1.yaxis.label.set_color('black')
# # ax1.tick_params(axis='y', which='both', colors='black')
# # ax1.set_ylim(0,)

# # # Plot the axis labels
# # ax1.set_ylabel('ALtitude (m)', fontsize=12)
# # ax1.set_xlabel('Date', fontsize=12, labelpad=-0.1)
# # plt.legend(bbox_to_anchor=(0.75, 0.95), loc=2, borderaxespad=0.)
