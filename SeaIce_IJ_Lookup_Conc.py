#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:23:56 2019

@author: ncp532
"""

# Date and Time handling package
from datetime import datetime	# functions to handle date and time

# FILE SYSTEM PACKAGES
from netCDF4 import Dataset				# function used to open multiple netcdf files

# DATA HANDLING PACKAGES
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------
# DEFINE THE DATASET

#-------------------------------
# SEA ICE SATELLITE DATA
 
# xarray
# 2018-19
#cubes2 = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_icdr_sh_f18_20180201_v01r00.nc',decode_cf=False,engine='netcdf4')
# 2017
#cubes2 = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20171231_v03r01.nc',decode_cf=False,engine='netcdf4')

# netcdf
# 2018-19
#cubes2 = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/seaice_conc_daily_icdr_sh_f18_20181120_v01r00.nc')
# 2017
#cubes2 = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20171231_v03r01.nc')

#-------------------------------
# IJ VALUES FOR THE LOOKUP FUNCTION

# CAMMPCAN 2017-18
#V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V1_17_M_ij.csv')
#V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V2_17_M_ij.csv')
#V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V3_17_M_ij.csv')
#
## CAMMPCAN 2018-19
#V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1_18_M_ij.csv')
#V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V2_18_M_ij.csv')
#V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V3_18_M_ij.csv')

# CAMMPCAN 2012
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/SIPEXII_M_ij.csv')

#------------------------------------------------------------------------------
# DEFINE THE NSIDC VARIABLES

#-------------------------------
# SEA ICE SATELLITE DATA

# NETCDF
#seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#y = cubes2.variables['latitude'][:] # y = 332, x = 316
#x = cubes2.variables['longitude'][:] # y = 332, x = 316

# XARRAY
#seaice_data = cubes2.seaice_conc_cdr[0,:,:] # y = 332, x = 316
#y = cubes2.latitude
#x = cubes2.longitude

#-------------------------------
# IJ VALUES FOR THE LOOKUP FUNCTION

## Date, Lat and Long
#Date_V1_17 = np.array(V1_17['Date_V1_17'])
#DateA = np.array(V1_17['Date'])
#iV1_17 = np.array(V1_17['I'])
#jV1_17 = np.array(V1_17['J'])
#
#Date_V2_17 = np.array(V2_17['Date_V2_17'])
#DateB = np.array(V2_17['Date'])
#iV2_17 = np.array(V2_17['I'])
#jV2_17 = np.array(V2_17['J'])
#
#Date_V3_17 = np.array(V3_17['Date_V3_17'])
#DateC = np.array(V3_17['Date'])
#iV3_17 = np.array(V3_17['I'])
#jV3_17 = np.array(V3_17['J'])
#
#Date_V1_18 = np.array(V1_18['Date_V1_18'])
#DateD = np.array(V1_18['Date'])
#iV1_18 = np.array(V1_18['I'])
#jV1_18 = np.array(V1_18['J'])
#
#Date_V2_18 = np.array(V2_18['Date_V2_18'])
#DateE = np.array(V2_18['Date'])
#iV2_18 = np.array(V2_18['I'])
#jV2_18 = np.array(V2_18['J'])
#
#Date_V3_18 = np.array(V3_18['Date_V3_18'])
#DateF = np.array(V3_18['Date'])
#iV3_18 = np.array(V3_18['I'])
#jV3_18 = np.array(V3_18['J'])

Date_SIPEXII = np.array(SIPEXII['Date_SIPEXII'])
DateG = np.array(SIPEXII['Date'])
iSIPEXII = np.array(SIPEXII['I'])
jSIPEXII = np.array(SIPEXII['J'])

#------------------------------------------------------------------------------
#CONVERT TO DATETIME FROM STRING

## V1_17
#date1=[]
#for i in range(len(Date_V1_17)):
#    date1.append(datetime.strptime(Date_V1_17[i],'%d/%m/%Y %H:%M:%S'))
#
## V2_17
#date2=[]
#for i in range(len(Date_V2_17)):
#    date2.append(datetime.strptime(Date_V2_17[i],'%d/%m/%Y %H:%M:%S'))
#
## V3_17
#date3=[]
#for i in range(len(Date_V3_17)):
#    date3.append(datetime.strptime(Date_V3_17[i],'%d/%m/%Y %H:%M:%S'))
#
## V1_18
#date4=[]
#for i in range(len(Date_V1_18)):
#    date4.append(datetime.strptime(Date_V1_18[i],'%d/%m/%Y %H:%M:%S'))
#
## V2_18
#date5=[]
#for i in range(len(Date_V2_18)):
#    date5.append(datetime.strptime(Date_V2_18[i],'%d/%m/%Y %H:%M:%S'))
#
## V3_18
#date7=[]
#for i in range(len(Date_V3_18)):
#    date7.append(datetime.strptime(Date_V3_18[i],'%d/%m/%Y %H:%M:%S'))

# SIPEXII
date8=[]
for i in range(len(Date_SIPEXII)):
    date8.append(datetime.strptime(Date_SIPEXII[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------------------------------------------------
# Set up reservoir arrays
#V1_17_SI = np.zeros(len(iV1_17))
#V2_17_SI = np.zeros(len(iV2_17))
#V3_17_SI = np.zeros(len(iV3_17))
#V1_18_SI = np.zeros(len(iV1_18))
#V2_18_SI = np.zeros(len(iV2_18))
#V3_18_SI = np.zeros(len(iV3_18))
SIPEXII_SI = np.zeros(len(iSIPEXII))

#------------------------------------------------------------------------------
# LOOK UP SEA ICE CONCENTRATION

##-------------------------------
## V1_17
#
## Set the folder location
#data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/'
#
## Loop over the dates in the voyage
#for i in range(len(DateA)):
#    # Set the file name
#    file = 'seaice_conc_daily_sh_f17_' + str(DateA[i]) + '_v03r01.nc'
#    # Select the file to open based on the date
#    cubes2 = Dataset(data_folder + file)
#    # Select the variable for sea ice concentration
#    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#    # Look up the sea ice concentration for the relevant grid square
#    V1_17_SI[i] = seaice_data[iV1_17[i]][jV1_17[i]] 
#
## Build a pandas dataframe
#dfV1_17_SI = np.column_stack((Date_V1_17,V1_17_SI))
#dfV1_17_SI = pd.DataFrame.from_dict(dfV1_17_SI)
#dfV1_17_SI.columns = ['Date_V1_17','Sea_Ice_Conc']
## Save the dataframe as a .csv
#dfV1_17_SI.to_csv('/Users/ncp532/Documents/Data/V1_17_M_SeaIce.csv')

##-------------------------------
## V2_17
#
## Set the folder location
#data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/'
#
## Loop over the dates in the voyage
#for i in range(len(DateB)):
#    # Set the file name
#    file = 'seaice_conc_daily_sh_f17_' + str(DateB[i]) + '_v03r01.nc'
#    # Select the file to open based on the date
#    cubes2 = Dataset(data_folder + file)
#    # Select the variable for sea ice concentration
#    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#    # Look up the sea ice concentration for the relevant grid square
#    V2_17_SI[i] = seaice_data[iV2_17[i]][jV2_17[i]] 
#
## Build a pandas dataframe
#dfV2_17_SI = np.column_stack((Date_V2_17,V2_17_SI))
#dfV2_17_SI = pd.DataFrame.from_dict(dfV2_17_SI)
#dfV2_17_SI.columns = ['Date_V2_17','Sea_Ice_Conc']
## Save the dataframe as a .csv
#dfV2_17_SI.to_csv('/Users/ncp532/Documents/Data/V2_17_M_SeaIce.csv')

##-------------------------------
## V3_17
#
## Set the folder location
#data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/'
#
## Loop over the dates in the voyage
#for i in range(len(DateC)):
#    # Set the file name
#    file = 'seaice_conc_daily_icdr_sh_f18_' + str(DateC[i]) + '_v01r00.nc'
#    # Select the file to open based on the date
#    cubes2 = Dataset(data_folder + file)
#    # Select the variable for sea ice concentration
#    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#    # Look up the sea ice concentration for the relevant grid square
#    V3_17_SI[i] = seaice_data[iV3_17[i]][jV3_17[i]] 
#
## Build a pandas dataframe
#dfV3_17_SI = np.column_stack((Date_V3_17,V3_17_SI))
#dfV3_17_SI = pd.DataFrame.from_dict(dfV3_17_SI)
#dfV3_17_SI.columns = ['Date_V3_17','Sea_Ice_Conc']
## Save the dataframe as a .csv
#dfV3_17_SI.to_csv('/Users/ncp532/Documents/Data/V3_17_M_SeaIce.csv')

##-------------------------------
## V1_18
#
## Set the folder location
#data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/'
#
## Loop over the dates in the voyage
#for i in range(len(DateD)):
#    # Set the file name
#    file = 'seaice_conc_daily_icdr_sh_f18_' + str(DateD[i]) + '_v01r00.nc'
#    # Select the file to open based on the date
#    cubes2 = Dataset(data_folder + file)
#    # Select the variable for sea ice concentration
#    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#    # Look up the sea ice concentration for the relevant grid square
#    V1_18_SI[i] = seaice_data[iV1_18[i]][jV1_18[i]] 
#
## Build a pandas dataframe
#dfV1_18_SI = np.column_stack((Date_V1_18,V1_18_SI))
#dfV1_18_SI = pd.DataFrame.from_dict(dfV1_18_SI)
#dfV1_18_SI.columns = ['Date_V1_18','Sea_Ice_Conc']
## Save the dataframe as a .csv
#dfV1_18_SI.to_csv('/Users/ncp532/Documents/Data/V1_18_M_SeaIce.csv')

##-------------------------------
## V2_18
#
## Set the folder location
#data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/'
#
## Loop over the dates in the voyage
#for i in range(len(DateE)):
#    # Set the file name
#    file = 'seaice_conc_daily_icdr_sh_f18_' + str(DateE[i]) + '_v01r00.nc'
#    # Select the file to open based on the date
#    cubes2 = Dataset(data_folder + file)
#    # Select the variable for sea ice concentration
#    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#    # Look up the sea ice concentration for the relevant grid square
#    V2_18_SI[i] = seaice_data[iV2_18[i]][jV2_18[i]] 
#
## Build a pandas dataframe
#dfV2_18_SI = np.column_stack((Date_V2_18,V2_18_SI))
#dfV2_18_SI = pd.DataFrame.from_dict(dfV2_18_SI)
#dfV2_18_SI.columns = ['Date_V2_18','Sea_Ice_Conc']
## Save the dataframe as a .csv
#dfV2_18_SI.to_csv('/Users/ncp532/Documents/Data/V2_18_M_SeaIce.csv')

##-------------------------------
## V3_18
#
## Set the folder location
#data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/'
#
## Loop over the dates in the voyage
#for i in range(len(DateF)):
#    # Set the file name
#    file = 'seaice_conc_daily_icdr_sh_f18_' + str(DateF[i]) + '_v01r00.nc'
#    # Select the file to open based on the date
#    cubes2 = Dataset(data_folder + file)
#    # Select the variable for sea ice concentration
#    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#    # Look up the sea ice concentration for the relevant grid square
#    V3_18_SI[i] = seaice_data[iV3_18[i]][jV3_18[i]] 
#
## Build a pandas dataframe
#dfV3_18_SI = np.column_stack((Date_V3_18,V3_18_SI))
#dfV3_18_SI = pd.DataFrame.from_dict(dfV3_18_SI)
#dfV3_18_SI.columns = ['Date_V3_18','Sea_Ice_Conc']
## Save the dataframe as a .csv
#dfV3_18_SI.to_csv('/Users/ncp532/Documents/Data/V3_18_M_SeaIce.csv')

#-------------------------------
# SIPEXII

# Set the folder location
data_folder = '/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2012/'

# Loop over the dates in the voyage
for i in range(len(DateG)):
    # Set the file name
    file = 'seaice_conc_daily_sh_f17_' + str(DateG[i]) + '_v03r01.nc'
    # Select the file to open based on the date
    cubes2 = Dataset(data_folder + file)
    # Select the variable for sea ice concentration
    seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2012
    # Look up the sea ice concentration for the relevant grid square
    SIPEXII_SI[i] = seaice_data[iSIPEXII[i]][jSIPEXII[i]] 

# Build a pandas dataframe
dfSIPEXII_SI = np.column_stack((Date_SIPEXII,SIPEXII_SI))
dfSIPEXII_SI = pd.DataFrame.from_dict(dfSIPEXII_SI)
dfSIPEXII_SI.columns = ['Date_SIPEXII','Sea_Ice_Conc']
# Save the dataframe as a .csv
dfSIPEXII_SI.to_csv('/Users/ncp532/Documents/Data/SIPEXII_M_SeaIce.csv')