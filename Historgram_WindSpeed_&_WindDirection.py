#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:27:01 2019

@author: ncp532
"""
# DESCRIPTION (VECTOR MEAN WIND DIRECTION)
# In terms of a description of the vector mean wind direction calculation
# the equation we use is:
 
# WD_vect = WD_s*WS_s/(WS_s+WS_p) + WD_p*WS_p/(WS_s+WS_p)
 
# Where:
# WD_vect is the vector mean wind direction
# WD_s and WD_p is the startboard and port wind direction respectively
# WS_s and WS_p is the starboard and port wind speed respectively.

# So basically, what we’re doing is weighting the sum of the wind directions
# from each sensor (port and starboard) by the wind speed measured by each
# sensor. Why its called a vector, is that the wind direction is the angle of
# the vector, while the wind speed is the vector’s magnitude.

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
from windrose import WindroseAxes

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd

#------------------------------------------------------------------------------
# Define the datasets

# CAMMPCAN (2017-18)
V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V01/CAMMPCAN_V1_underway_60.csv') # V1 Cruise data
V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V02/CAMMPCAN_V2_underway_60.csv') # V2 Cruise data 
V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V03/CAMMPCAN_V3_underway_60.csv') # V3 Cruise data
V4_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V04/CAMMPCAN_V4_underway_60.csv') # V4 Cruise data

# CAMMPCAN (2018-19)
V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V01/CAMMPCAN_V1_underway_60.csv') # V1 Cruise data
V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V02/CAMMPCAN_V2_underway_60.csv') # V2 Cruise data 
V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V03/CAMMPCAN_V3_underway_60.csv') # V3 Cruise data
V4_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V04/CAMMPCAN_V4_underway_60.csv') # V4 Cruise data

# PCAN (2017)
PCAN = pd.read_csv('/Users/ncp532/Documents/Data/PCAN_2017/PCAN_underway_data.csv') # PCAN Cruise data

# SIPEXII (2012)
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII/SIPEXII_underway_data.csv') # SIPEXII Cruise data

#------------------------------------------------------------------------------
# FILTER OUT DATA FROM OVER 60 DEGREES SOUTH

# CAMMPCAN (2017-18)
LATS = V1_17['LATITUDE']<=-60
V1_17_60S = V1_17[LATS]

LATS = V2_17['LATITUDE']<=-60
V2_17_60S = V2_17[LATS]

LATS = V3_17['LATITUDE']<=-60
V3_17_60S = V3_17[LATS]

# CAMMPCAN (2018-19)
LATS = V1_18['LATITUDE']<=-60
V1_18_60S = V1_18[LATS]

LATS = V2_18['LATITUDE']<=-60
V2_18_60S = V2_18[LATS]

LATS = V3_18['LATITUDE']<=-60
V3_18_60S = V3_18[LATS]

# PCAN
LATS = PCAN['Latitude (degree_north)']<=-60
PCAN_60S = PCAN[LATS]

# SIPEXII
LATS = SIPEXII['LATITUDE_DEGNORTH']<=-60
SIPEXII_60S = SIPEXII[LATS]

#------------------------------------------------------------------------------
# Define the variables

# CAMMPCAN (2017-18)
# V1----------------
# Wind Direction
WD_s1_17 = np.array(V1_17['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p1_17 = np.array(V1_17['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s1_17 = np.array(V1_17['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s1_17 = WS_s1_17 * 0.514444444 # Convert from knots to m/s
WS_p1_17 = np.array(V1_17['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p1_17 = WS_p1_17 * 0.514444444 # Convert from knots to m/s
WS_1_17 = (WS_s1_17 + WS_p1_17)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect1_17 = ((WD_s1_17 * WS_s1_17) / (WS_s1_17 + WS_p1_17)) + ((WD_p1_17 * WS_p1_17) / (WS_s1_17 + WS_p1_17)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df1_17 = np.column_stack((WD_vect1_17,WS_1_17))
df1_17 = pd.DataFrame.from_dict(df1_17)
df1_17.columns = ['WD_vect1_17', 'WS_1_17']

# V2----------------
# Wind Direction
WD_s2_17 = np.array(V2_17['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p2_17 = np.array(V2_17['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s2_17 = np.array(V2_17['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s2_17 = WS_s2_17 * 0.514444444 # Convert from knots to m/s
WS_p2_17 = np.array(V2_17['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p2_17 = WS_p2_17 * 0.514444444 # Convert from knots to m/s
WS_2_17 = (WS_s2_17 + WS_p2_17)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect2_17 = ((WD_s2_17 * WS_s2_17) / (WS_s2_17 + WS_p2_17)) + ((WD_p2_17 * WS_p2_17) / (WS_s2_17 + WS_p2_17)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df2_17 = np.column_stack((WD_vect2_17,WS_2_17))
df2_17 = pd.DataFrame.from_dict(df2_17)
df2_17.columns = ['WD_vect2_17', 'WS_2_17']

# V3----------------
# Wind Direction
WD_s3_17 = np.array(V3_17['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p3_17 = np.array(V3_17['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s3_17 = np.array(V3_17['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s3_17 = WS_s3_17 * 0.514444444 # Convert from knots to m/s
WS_p3_17 = np.array(V3_17['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p3_17 = WS_p3_17 * 0.514444444 # Convert from knots to m/s
WS_3_17 = (WS_s3_17 + WS_p3_17)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect3_17 = ((WD_s3_17 * WS_s3_17) / (WS_s3_17 + WS_p3_17)) + ((WD_p3_17 * WS_p3_17) / (WS_s3_17 + WS_p3_17)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df3_17 = np.column_stack((WD_vect3_17,WS_3_17))
df3_17 = pd.DataFrame.from_dict(df3_17)
df3_17.columns = ['WD_vect3_17', 'WS_3_17']

# V4----------------
# Wind Direction
WD_s4_17 = np.array(V4_17['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p4_17 = np.array(V4_17['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s4_17 = np.array(V4_17['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s4_17 = WS_s4_17 * 0.514444444 # Convert from knots to m/s
WS_p4_17 = np.array(V4_17['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p4_17 = WS_p4_17 * 0.514444444 # Convert from knots to m/s
WS_4_17 = (WS_s4_17 + WS_p4_17)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect4_17 = ((WD_s4_17 * WS_s4_17) / (WS_s4_17 + WS_p4_17)) + ((WD_p4_17 * WS_p4_17) / (WS_s4_17 + WS_p4_17)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df4_17 = np.column_stack((WD_vect4_17,WS_4_17))
df4_17 = pd.DataFrame.from_dict(df4_17)
df4_17.columns = ['WD_vect4_17', 'WS_4_17']

# CAMMPCAN (2018-19)
# V1----------------
# Wind Direction
WD_s1_18 = np.array(V1_18['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p1_18 = np.array(V1_18['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s1_18 = np.array(V1_18['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s1_18 = WS_s1_18 * 0.514444444 # Convert from knots to m/s
WS_p1_18 = np.array(V1_18['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p1_18 = WS_p1_18 * 0.514444444 # Convert from knots to m/s
WS_1_18 = (WS_s1_18 + WS_p1_18)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect1_18 = ((WD_s1_18 * WS_s1_18) / (WS_s1_18 + WS_p1_18)) + ((WD_p1_18 * WS_p1_18) / (WS_s1_18 + WS_p1_18)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df1_18 = np.column_stack((WD_vect1_18,WS_1_18))
df1_18 = pd.DataFrame.from_dict(df1_18)
df1_18.columns = ['WD_vect1_18', 'WS_1_18']

# V2----------------
# Wind Direction
WD_s2_18 = np.array(V2_18['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p2_18 = np.array(V2_18['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s2_18 = np.array(V2_18['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s2_18 = WS_s2_18 * 0.514444444 # Convert from knots to m/s
WS_p2_18 = np.array(V2_18['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p2_18 = WS_p2_18 * 0.514444444 # Convert from knots to m/s
WS_2_18 = (WS_s2_18 + WS_p2_18)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect2_18 = ((WD_s2_18 * WS_s2_18) / (WS_s2_18 + WS_p2_18)) + ((WD_p2_18 * WS_p2_18) / (WS_s2_18 + WS_p2_18)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df2_18 = np.column_stack((WD_vect2_18,WS_2_18))
df2_18 = pd.DataFrame.from_dict(df2_18)
df2_18.columns = ['WD_vect2_18', 'WS_2_18']

# V3----------------
# Wind Direction
WD_s3_18 = np.array(V3_18['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p3_18 = np.array(V3_18['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s3_18 = np.array(V3_18['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s3_18 = WS_s3_18 * 0.514444444 # Convert from knots to m/s
WS_p3_18 = np.array(V3_18['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p3_18 = WS_p3_18 * 0.514444444 # Convert from knots to m/s
WS_3_18 = (WS_s3_18 + WS_p3_18)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect3_18 = ((WD_s3_18 * WS_s3_18) / (WS_s3_18 + WS_p3_18)) + ((WD_p3_18 * WS_p3_18) / (WS_s3_18 + WS_p3_18)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df3_18 = np.column_stack((WD_vect3_18,WS_3_18))
df3_18 = pd.DataFrame.from_dict(df3_18)
df3_18.columns = ['WD_vect3_18', 'WS_3_18']

# V4----------------
# Wind Direction
WD_s4_18 = np.array(V4_18['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_p4_18 = np.array(V4_18['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_s4_18 = np.array(V4_18['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_s4_18 = WS_s4_18 * 0.514444444 # Convert from knots to m/s
WS_p4_18 = np.array(V4_18['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_p4_18 = WS_p4_18 * 0.514444444 # Convert from knots to m/s
WS_4_18 = (WS_s4_18 + WS_p4_18)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect4_18 = ((WD_s4_18 * WS_s4_18) / (WS_s4_18 + WS_p4_18)) + ((WD_p4_18 * WS_p4_18) / (WS_s4_18 + WS_p4_18)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
df4_18 = np.column_stack((WD_vect4_18,WS_4_18))
df4_18 = pd.DataFrame.from_dict(df4_18)
df4_18.columns = ['WD_vect4_18', 'WS_4_18']

# PCAN (2017)--------
# Wind Direction
WD_sPCAN = np.array(PCAN['Starboard True Wind Direction (degree)']) # starboard side wind direction (correlated)
WD_pPCAN = np.array(PCAN['Port True Wind Direction (degree)']) # port side wind direction (correlated)
# Wind Speed
WS_sPCAN = np.array(PCAN['Starboard True Wind Speed (knot)']) # starboard side wind speed (correlated)
WS_sPCAN = WS_sPCAN * 0.514444444 # Convert from knots to m/s
WS_pPCAN = np.array(PCAN['Port True Wind Speed (knot)']) # port side wind speed (correlated)
WS_pPCAN = WS_pPCAN * 0.514444444 # Convert from knots to m/s
WS_PCAN = (WS_sPCAN + WS_pPCAN)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vectPCAN = ((WD_sPCAN * WS_sPCAN) / (WS_sPCAN + WS_pPCAN)) + ((WD_pPCAN * WS_pPCAN) / (WS_sPCAN + WS_pPCAN)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
dfPCAN = np.column_stack((WD_vectPCAN,WS_PCAN))
dfPCAN = pd.DataFrame.from_dict(dfPCAN)
dfPCAN.columns = ['WD_vectPCAN', 'WS_PCAN']

# SIPEXII (2012)----
# Wind Direction
WD_sSIP = np.array(SIPEXII['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
WD_pSIP = np.array(SIPEXII['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
# Wind Speed
WS_sSIP = np.array(SIPEXII['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
WS_sSIP = WS_sSIP * 0.514444444 # Convert from knots to m/s
WS_pSIP = np.array(SIPEXII['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
WS_pSIP = WS_pSIP * 0.514444444 # Convert from knots to m/s
WS_SIP = (WS_sSIP + WS_pSIP)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vectSIP = ((WD_sSIP * WS_sSIP) / (WS_sSIP + WS_pSIP)) + ((WD_pSIP * WS_pSIP) / (WS_sSIP + WS_pSIP)) # Calculate the vector mean wind direction
# Build DataFrame of WindDirection & WindSpeed
dfSIPII = np.column_stack((WD_vectSIP,WS_SIP))
dfSIPII = pd.DataFrame.from_dict(dfSIPII)
dfSIPII.columns = ['WD_vectSIP', 'WS_SIP']
#------------------------------------------------------------------------------
# PLOT HISTOGRAMS OF WIND SPEED AND VECTOR MEAN WIND DIRECTION

#---------------------------------
# CAMMPCAN (2017-18)
fig1 = plt.figure(figsize=(20,15)) # (width, height)

rect_1 = [0.04, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig1.add_axes(rect_1)    # creates the axes of specified dimensions
df1_17['WD_vect1_17'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2017-18 V1", position=(0.5, 1.1))

rect_2 = [0.04, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig1.add_axes(rect_2)    # creates the axes of specified dimensions
df1_17['WS_1_17'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_3 = [-0.01, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig1, rect_3)    # creates the axes of specified dimensions 
fig1.add_axes(ax)
ax.bar(WD_vect1_17, WS_1_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_4 = [0.29, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig1.add_axes(rect_4)    # creates the axes of specified dimensions
df2_17['WD_vect2_17'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2017-18 V2", position=(0.5, 1.1))

rect_5 = [0.29, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig1.add_axes(rect_5)    # creates the axes of specified dimensions
df2_17['WS_2_17'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_6 = [0.24, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig1, rect_6)
fig1.add_axes(ax)
ax.bar(WD_vect2_17, WS_2_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_7 = [0.54, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig1.add_axes(rect_7)    # creates the axes of specified dimensions
df3_17['WD_vect3_17'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2017-18 V3", position=(0.5, 1.1))

rect_8 = [0.54, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig1.add_axes(rect_8)    # creates the axes of specified dimensions
df3_17['WS_3_17'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_9 = [0.49, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig1, rect_9)
fig1.add_axes(ax)
ax.bar(WD_vect3_17, WS_3_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_10 = [0.79, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig1.add_axes(rect_10)    # creates the axes of specified dimensions
df4_17['WD_vect4_17'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2017-18 V4", position=(0.5, 1.1))

rect_11 = [0.79, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig1.add_axes(rect_11)    # creates the axes of specified dimensions
df4_17['WS_4_17'].hist(bins=100, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlim(0,27)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_12 = [0.74, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig1, rect_12)
fig1.add_axes(ax)
ax.bar(WD_vect4_17, WS_4_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

#---------------------------------
# CAMMPCAN (2018-19)
fig2 = plt.figure(figsize=(20,15)) # (width, height)

rect_1 = [0.04, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig2.add_axes(rect_1)    # creates the axes of specified dimensions
df1_18['WD_vect1_18'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2018-19 V1", position=(0.5, 1.1))

rect_2 = [0.04, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig2.add_axes(rect_2)    # creates the axes of specified dimensions
df1_18['WS_1_18'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_3 = [-0.01, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig2, rect_3)    # creates the axes of specified dimensions 
fig2.add_axes(ax)
ax.bar(WD_vect1_18, WS_1_18, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_4 = [0.29, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig2.add_axes(rect_4)    # creates the axes of specified dimensions
df2_18['WD_vect2_18'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2018-19 V2", position=(0.5, 1.1))

rect_5 = [0.29, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig2.add_axes(rect_5)    # creates the axes of specified dimensions
df2_18['WS_2_18'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_6 = [0.24, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig2, rect_6)
fig2.add_axes(ax)
ax.bar(WD_vect2_18, WS_2_18, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_7 = [0.54, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig2.add_axes(rect_7)    # creates the axes of specified dimensions
df3_18['WD_vect3_18'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2018-19 V3", position=(0.5, 1.1))

rect_8 = [0.54, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig2.add_axes(rect_8)    # creates the axes of specified dimensions
df3_18['WS_3_18'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_9 = [0.49, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig2, rect_9)
fig2.add_axes(ax)
ax.bar(WD_vect3_18, WS_3_18, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_10 = [0.79, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig2.add_axes(rect_10)    # creates the axes of specified dimensions
df4_18['WD_vect4_18'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("CAMMPCAN 2018-19 V4", position=(0.5, 1.1))

rect_11 = [0.79, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig2.add_axes(rect_11)    # creates the axes of specified dimensions
df4_18['WS_4_18'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_12 = [0.74, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig2, rect_12)
fig2.add_axes(ax)
ax.bar(WD_vect4_18, WS_4_18, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

#---------------------------------
# PCAN & SIPEXII
fig3 = plt.figure(figsize=(20,15)) # (width, height)

rect_1 = [0.29, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig3.add_axes(rect_1)    # creates the axes of specified dimensions
dfPCAN['WD_vectPCAN'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("PCAN (2017)", position=(0.5, 1.1))

rect_2 = [0.29, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig3.add_axes(rect_2)    # creates the axes of specified dimensions
dfPCAN['WS_PCAN'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_3 = [0.24, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig3, rect_3)    # creates the axes of specified dimensions 
fig3.add_axes(ax)
ax.bar(WD_vectPCAN, WS_PCAN, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))

rect_4 = [0.54, 0.7, 0.2, 0.2]     # [left, bottom, width, height] as a fraction of total figure size
ax = fig3.add_axes(rect_4)    # creates the axes of specified dimensions
dfSIPII['WD_vectSIP'].hist(bins=10, edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind direction")
ax.set_ylabel("Frequency")
ax.set_title("SIPEXII (2012)", position=(0.5, 1.1))

rect_5 = [0.54, 0.375, 0.2, 0.2]     # [left, bottom, width, height]
ax = fig3.add_axes(rect_5)    # creates the axes of specified dimensions
dfSIPII['WS_SIP'].hist(bins=10, color = "orange", edgecolor = 'black', linewidth = 0.5)
ax.set_xlabel("wind speed (m/s)")
ax.set_ylabel("Frequency")

rect_6 = [0.49, 0.05, 0.225, 0.225]
ax = WindroseAxes(fig3, rect_6)
fig3.add_axes(ax)
ax.bar(WD_vectSIP, WS_SIP, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
ax.legend(title="wind speed (m/s)", loc=(1.1, 0))