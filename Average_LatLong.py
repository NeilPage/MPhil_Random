#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 11:20:00 2019

@author: ncp532
"""

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
#from windrose import WindroseAxes

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
import itertools

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

#------------------------------------------------------------------------------
# OBSERVATIONAL DATA

# CAMMPCAN 2017-18
V1_17_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V01/CAMMPCAN_V1_underway_60.csv') 
V2_17_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V02/CAMMPCAN_V2_underway_60.csv') 
V3_17_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V03/CAMMPCAN_V3_underway_60.csv') 

# CAMMPCAN 2018-19
V1_18_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V01/CAMMPCAN_V1_underway_60.csv') 
V2_18_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V02/CAMMPCAN_V2_underway_60.csv') 
V3_18_Met = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V03/CAMMPCAN_V3_underway_60.csv') 

# SIPEXII 2012
SIPEXII_Met = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/201213001.csv') #SIPEXII_underway_60.csv')

#------------------------------------------------------------------------------
# Define the variables

# CAMMPCAN (2017-18)
Lat_V1_17 = np.array(V1_17_Met['LATITUDE'])
Lon_V1_17 = np.array(V1_17_Met['LONGITUDE'])

Lat_V2_17 = np.array(V2_17_Met['LATITUDE'])
Lon_V2_17 = np.array(V2_17_Met['LONGITUDE'])

Lat_V3_17 = np.array(V3_17_Met['LATITUDE'])
Lon_V3_17 = np.array(V3_17_Met['LONGITUDE'])

# CAMMPCAN (2018-19)
Lat_V1_18 = np.array(V1_18_Met['LATITUDE'])
Lon_V1_18 = np.array(V1_18_Met['LONGITUDE'])

Lat_V2_18 = np.array(V2_18_Met['LATITUDE'])
Lon_V2_18 = np.array(V2_18_Met['LONGITUDE'])

Lat_V3_18 = np.array(V3_18_Met['LATITUDE'])
Lon_V3_18 = np.array(V3_18_Met['LONGITUDE'])

# SIPEXII (2012)
Lat_SIPEXII = np.array(SIPEXII_Met['LATITUDE'])
Lon_SIPEXII = np.array(SIPEXII_Met['LONGITUDE'])

#------------------------------------------------------------------------------
# SET THE DATE AND TIME
#------------------------------------
# V1_17
datM1 = np.array(V1_17_Met['Date'])
timM1 = np.array(V1_17_Met['Time'])
dattimM1 = datM1+' '+timM1

#CONVERT TO DATETIME FROM STRING
dateM1=[]
for i in range(len(dattimM1)):
    dateM1.append(datetime.strptime(dattimM1[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------    
# V2_17
datM2 = np.array(V2_17_Met['Date'])
timM2 = np.array(V2_17_Met['Time'])
dattimM2 = datM2+' '+timM2

#CONVERT TO DATETIME FROM STRING
dateM2=[]
for i in range(len(dattimM2)):
    dateM2.append(datetime.strptime(dattimM2[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------
# V3_17
datM3 = np.array(V3_17_Met['Date'])
timM3 = np.array(V3_17_Met['Time'])
dattimM3 = datM3+' '+timM3

#CONVERT TO DATETIME FROM STRING
dateM3=[]
for i in range(len(dattimM3)):
    dateM3.append(datetime.strptime(dattimM3[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------
# V1_18
datM4 = np.array(V1_18_Met['Date'])
timM4 = np.array(V1_18_Met['Time'])
dattimM4 = datM4+' '+timM4

#CONVERT TO DATETIME FROM STRING
dateM4=[]
for i in range(len(dattimM4)):
    dateM4.append(datetime.strptime(dattimM4[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------
# V2_18
datM5 = np.array(V2_18_Met['Date'])
timM5 = np.array(V2_18_Met['Time'])
dattimM5 = datM5+' '+timM5

#CONVERT TO DATETIME FROM STRING
dateM5=[]
for i in range(len(dattimM5)):
    dateM5.append(datetime.strptime(dattimM5[i],'%d/%m/%Y %H:%M:%S'))
    
#------------------------------------
# V3_18
datM6 = np.array(V3_18_Met['Date'])
timM6 = np.array(V3_18_Met['Time'])
dattimM6 = datM6+' '+timM6

#CONVERT TO DATETIME FROM STRING
dateM6=[]
for i in range(len(dattimM6)):
    dateM6.append(datetime.strptime(dattimM6[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------
# SIPEXII
datM7 = np.array(SIPEXII_Met['Date'])
timM7 = np.array(SIPEXII_Met['Time'])
dattimM7 = datM7+' '+timM7

#CONVERT TO DATETIME FROM STRING
dateM7=[]
for i in range(len(dattimM7)):
    dateM7.append(datetime.strptime(dattimM7[i],'%d/%m/%Y %H:%M:%S'))

#------------------------------------------------------------------------------
# CALCULATE THE BRO DAILY AVERAGE (MEAN)

# Function to calculate the daily mean
def dailyM(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('T').mean()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 

# BrO Daily Means

# Latitude Daily Means
Lat_V1_17_DM, dateM1_DM=dailyM(Lat_V1_17[:],dateM1) # V1_17
Lat_V2_17_DM, dateM2_DM=dailyM(Lat_V2_17[:],dateM2) # V2_17
Lat_V3_17_DM, dateM3_DM=dailyM(Lat_V3_17[:],dateM3) # V3_17
Lat_V1_18_DM, dateM4_DM=dailyM(Lat_V1_18[:],dateM4) # V1_18
Lat_V2_18_DM, dateM5_DM=dailyM(Lat_V2_18[:],dateM5) # V2_18
Lat_V3_18_DM, dateM6_DM=dailyM(Lat_V3_18[:],dateM6) # V3_18
Lat_SIPEXII_DM, dateM7_DM=dailyM(Lat_SIPEXII[:],dateM7) # V3_18

# Longitude Daily Means
Lon_V1_17_DM, dateM1_DM=dailyM(Lon_V1_17[:],dateM1) # V1_17
Lon_V2_17_DM, dateM2_DM=dailyM(Lon_V2_17[:],dateM2) # V2_17
Lon_V3_17_DM, dateM3_DM=dailyM(Lon_V3_17[:],dateM3) # V3_17
Lon_V1_18_DM, dateM4_DM=dailyM(Lon_V1_18[:],dateM4) # V1_18
Lon_V2_18_DM, dateM5_DM=dailyM(Lon_V2_18[:],dateM5) # V2_18
Lon_V3_18_DM, dateM6_DM=dailyM(Lon_V3_18[:],dateM6) # V3_18
Lon_SIPEXII_DM, dateM7_DM=dailyM(Lon_SIPEXII[:],dateM7) # V3_18

#------------------------------------------------------------------------------
# BUILD A NEW DATAFRAME FOR LATITUDE AND LONGITUDE DAILY AVERAGE AND SAVE AS A .CSV 

# Hourly/Daily
#LatLon = [dateM1_DM,Lat_V1_17_DM,Lon_V1_17_DM,dateM2_DM,Lat_V2_17_DM,Lon_V2_17_DM,dateM3_DM,Lat_V3_17_DM,Lon_V3_17_DM,dateM4_DM,Lat_V1_18_DM,Lon_V1_18_DM,dateM5_DM,Lat_V2_18_DM,Lon_V2_18_DM,dateM6_DM,Lat_V3_18_DM,Lon_V3_18_DM]
#Minute
LatLon = [dateM1,Lat_V1_17,Lon_V1_17,dateM2,Lat_V2_17,Lon_V2_17,dateM3,Lat_V3_17,Lon_V3_17,dateM4,Lat_V1_18,Lon_V1_18,dateM5,Lat_V2_18,Lon_V2_18,dateM6,Lat_V3_18,Lon_V3_18,dateM7_DM,Lat_SIPEXII_DM,Lon_SIPEXII_DM]


# Function to format arrays to the same length
def stack_padding(l):
    return np.column_stack((itertools.zip_longest(*l, fillvalue=0)))

# Apply the function and transpose
dfLatLon = stack_padding(LatLon)
dfLatLon = np.transpose(dfLatLon)

# Convert to a Pandas dataframe and add column names
dfLatLon = pd.DataFrame.from_dict(dfLatLon)
dfLatLon.columns = ['Date_V1_17','Lat_V1_17','Lon_V1_17','Date_V2_17','Lat_V2_17','Lon_V2_17','Date_V3_17','Lat_V3_17','Lon_V3_17','Date_V1_18','Lat_V1_18','Lon_V1_18','Date_V2_18','Lat_V2_18','Lon_V2_18','Date_V3_18','Lat_V3_18','Lon_V3_18','DateSIPEXII','Lat_SIPEXII','Lon_SIPEXII']

# Save as a .csv file
dfLatLon.to_csv('/Users/ncp532/Documents/Data/SIPEXII_Minute_LatLong.csv')
