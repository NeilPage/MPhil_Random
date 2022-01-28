#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:16:43 2019

@author: ncp532
"""
# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files
import xarray as xr                 # xarray to open multiple netcdf files    

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import interpolate

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# Exhaust Filter
exhaust=xr.open_dataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/Exhaust_ID/AAS_4292_ExhaustID_201718_AA_MARCUS.nc')

# NO2 Data
V1_NO2 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/TG_Retrieval_V1/V1_NO2_QAQC.csv')

# BrO Data
V1_BrO = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/TG_Retrieval_V1/V1_BrO_QAQC.csv')

# O3 & CO Data
V1_O3 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/ARM_O3_CO_N2O.csv')

# Met data
V1_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V01/CAMMPCAN_V1_underway_60.csv') 

# Hg Data
V1_Hg = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V1_Hg0_Lat_Long_17-18.csv')

#------------------------------------------------------------------------------
# Define the variables

# Latitude / Longitude
Lat_V1_18 = np.array(V1_Met['LATITUDE'])
Long_V1_18 = np.array(V1_Met['LONGITUDE'])

# Wind Direction
WD_s1_18 = np.array(V1_Met['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p1_18 = np.array(V1_Met['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)

# Wind Speed
WS_s1_18 = np.array(V1_Met['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s1_18 = WS_s1_18 * 0.514444444 # Convert from knots to m/s
WS_p1_18 = np.array(V1_Met['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p1_18 = WS_p1_18 * 0.514444444 # Convert from knots to m/s
WS_1_18 = (WS_s1_18 + WS_p1_18)/2 # Average the wind speed for port and starboard

# Vector Mean Wind Direction
WD_vect1_18 = ((WD_s1_18 * WS_s1_18) / (WS_s1_18 + WS_p1_18)) + ((WD_p1_18 * WS_p1_18) / (WS_s1_18 + WS_p1_18)) # Calculate the vector mean wind direction

# BrO surface volume mixing ratio (VMR)
BrO_vmr_V1_18 = np.array(V1_BrO['surf_vmr(ppmv)'])
BrO_vmr_V1_18 = BrO_vmr_V1_18 * 1e6 # convert from ppmv to pptv

# NO2 surface volume mixing ratio (VMR)
NO2_vmr_V1_18 = np.array(V1_NO2['surf_vmr(ppmv)'])
NO2_vmr_V1_18 = NO2_vmr_V1_18 * 1e3 # convert from ppmv to ppbv

# Solar Radiation (W/m2)
Sol_s1_18 = np.array(V1_Met['RAD_SLR_STRBRD_WPERM2']) # starboard side solar radiation
Sol_p1_18 = np.array(V1_Met['RAD_SLR_PORT_WPERM2']) # port side solar radiation
Sol_V1_18 = np.array(V1_Met['Rad_Average'])

# Temperature (C)
Temp_s1_18 = np.array(V1_Met['TEMP_AIR_STRBRD_DEGC']) # starboard side temperature
Temp_p1_18 = np.array(V1_Met['TEMP_AIR_PORT_DEGC']) # port side temperature
Temp_V1_18 = (Temp_s1_18 + Temp_p1_18)/2 # Average the temperature for port and starboard

# Relative Humidity (RH)
RH_s1_18 = np.array(V1_Met['REL_HUMIDITY_STRBRD_PERCENT']) # starboard side relative humidity (%)
RH_p1_18 = np.array(V1_Met['REL_HUMIDITY_PORT_PERCENT']) # port side relative humidity (%)
RH_V1_18 = (RH_s1_18 + RH_p1_18)/2 # Average the RH for port and starboard (%)

# Atmospheric pressure (hPa)
Pres_s1_18 = np.array(V1_Met['ATM_PRESS_HPA']) # Atmospheric Pressure (hPa)

# O3 (ppb)
O3_V1_18 = np.array(V1_O3['O3_(ppb)']) # all O3 data

# CO (ppb)
CO_V1_18 = np.array(V1_O3['CO_Dry']) # CO ()

# Hg Data
Hg0_V1_18 = np.array(V1_Hg['ng/m3']) # Hg0
RM_V1_18 = np.array(V1_Hg['RM_pg/m3']) # RM
RM_StDev = np.array(V1_Hg['RM_StDev']) # RM StDev

# MAD (Median Absolute Distribution)
M3_T0 = np.array(exhaust.exhaust_3mad00thresh) # Exhaust ID using 3 MAD and 0% window threshold.
M3_T5 = np.array(exhaust.exhaust_3mad05thresh) # Exhaust ID using 3 MAD and 5% window threshold.
M3_T10 = np.array(exhaust.exhaust_3mad10thresh) # Exhaust ID using 3 MAD and 10% window threshold.
M4_T0 = np.array(exhaust.exhaust_4mad00thresh) # Exhaust ID using 4 MAD and 0% window threshold.
M4_T1 = np.array(exhaust.exhaust_4mad01thresh) # Exhaust ID using 4 MAD and 1% window threshold.
M4_T2 = np.array(exhaust.exhaust_4mad02thresh) # Exhaust ID using 4 MAD and 2% window threshold.
M4_T3 = np.array(exhaust.exhaust_4mad03thresh) # Exhaust ID using 4 MAD and 3% window threshold.
M4_T4 = np.array(exhaust.exhaust_4mad04thresh) # Exhaust ID using 4 MAD and 4% window threshold.
M4_T5 = np.array(exhaust.exhaust_4mad05thresh) # Exhaust ID using 4 MAD and 5% window threshold.
M4_T10 = np.array(exhaust.exhaust_4mad10thresh) # Exhaust ID using 4 MAD and 10% window threshold.

#------------------------------------------------------------------------------
# SET THE DATE AND TIME

# V1_BrO
dat = np.array(V1_BrO['Date'])
tim = np.array(V1_BrO['Time'])
dattim = dat+' '+tim

#CONVERT TO DATETIME FROM STRING
date=[]
for i in range(len(dattim)):
    date.append(datetime.strptime(dattim[i],'%d/%m/%Y %H:%M:%S'))

# Insert DATETIME to the df
V1_BrO.insert(0, 'DateTime', date)

#----------------------
# V1_O3
dat2 = np.array(V1_O3['Date'])
tim2 = np.array(V1_O3['Time'])
dattim2 = dat2+' '+tim2

#CONVERT TO DATETIME FROM STRING
date2=[]
for i in range(len(dattim2)):
    date2.append(datetime.strptime(dattim2[i],'%d/%m/%Y %H:%M:%S'))

# Insert DATETIME to the df
V1_O3.insert(0, 'DateTime', date2)

#----------------------
# V1_Met
dat3 = np.array(V1_Met['Date'])
tim3 = np.array(V1_Met['Time'])
dattim3 = dat3+' '+tim3

#CONVERT TO DATETIME FROM STRING
date3=[]
for i in range(len(dattim3)):
    date3.append(datetime.strptime(dattim3[i],'%d/%m/%Y %H:%M:%S'))

# Insert DATETIME to the df
V1_Met.insert(0, 'DateTime', date3)

#----------------------
# V1_NO2
dat4 = np.array(V1_NO2['Date'])
tim4 = np.array(V1_NO2['Time'])
dattim4 = dat4+' '+tim4

#CONVERT TO DATETIME FROM STRING
date4=[]
for i in range(len(dattim4)):
    date4.append(datetime.strptime(dattim4[i],'%d/%m/%Y %H:%M:%S'))

# Insert DATETIME to the df
V1_NO2.insert(0, 'DateTime', date4)

#----------------------
# V1_Hg0
dat5 = np.array(V1_Hg['Date'])
tim5 = np.array(V1_Hg['Time'])
dattim5 = dat5+' '+tim5

#CONVERT TO DATETIME FROM STRING
date5=[]
for i in range(len(dattim5)):
    date5.append(datetime.strptime(dattim5[i],'%d/%m/%Y %H:%M:%S'))

# Insert DATETIME to the df
V1_Hg.insert(0, 'DateTime', date5)

#----------------------
# Exhaust Filter
Time6 = np.array(exhaust.time)
ts6 = pd.to_datetime((Time6))
ts6 = ts6.strftime('%d/%m/%Y %H:%M:%S')

#CONVERT TO DATETIME FROM STRING
date6=[]
for i in range(len(ts6)):
    date6.append(datetime.strptime(ts6[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------------------------------------------------
# Interpolate data from 1 min to 1 second

def second(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('1S').bfill()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 

# BrO (1 second resample)  
BrO_1sec, date_1sec=second(BrO_vmr_V1_18,date) # O3

df = V1_BrO
df['DT'] = pd.to_datetime(df['DateTime'])
df = df.set_index('DT')
df = df.resample('1S').bfill()

# BUILD A NEW DATAFRAME FOR IMPORT TO R AND SAVE AS A .CSV 

exh = np.column_stack((ts,M3_T0,M3_T5,M3_T10,M4_T0,M4_T1,M4_T2,M4_T3,M4_T4,M4_T5,M4_T10))
df1 = pd.DataFrame.from_dict(exh)
df1.columns = ['Time','M3_T0','M3_T5','M3_T10','M4_T0','M4_T1','M4_T2','M4_T3','M4_T4','M4_T5','M4_T10']

# Export the combined .csv
df1.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/DF1_exhaust.csv')