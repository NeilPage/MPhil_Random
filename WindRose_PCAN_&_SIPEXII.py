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

# PCAN (2017)
#PCAN = pd.read_csv('/Users/ncp532/Documents/Data/PCAN_2017/PCAN_underway_data.csv') # PCAN Cruise data

# SIPEXII (2012)
#SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/SIPEXII_underway_60.csv') # SIPEXII Cruise data
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/ShipTrack/SIPEXII_underway_60.csv', index_col=0) 

#------------------------------------------------------------------------------
# FILTER OUT DATA CORRESPONDING TO BELOW 65 DEGREES SOUTH
## PCAN
#LATS = PCAN['Latitude (degree_north)']<=-65
#PCAN_60S = PCAN[LATS]

# SIPEXII
LATS = SIPEXII['latitude']<=-65
SIPEXII_60S = SIPEXII[LATS]

## FILTER OUT DATA CORRESPONDING TO ABOVE 65 DEGREES SOUTH
## PCAN
#LATS = PCAN['Latitude (degree_north)']>-65
#PCAN_N60S = PCAN[LATS]

# SIPEXII
LATS = SIPEXII['latitude']>-64
SIPEXII_N60S = SIPEXII[LATS]

#------------------------------------------------------------------------------
# Define the variables

# ALL LATITUDES
#
## PCAN (2017)--------
## Wind Direction
#WD_sPCAN = np.array(PCAN['Starboard True Wind Direction (degree)']) # starboard side wind direction (correlated)
#WD_pPCAN = np.array(PCAN['Port True Wind Direction (degree)']) # port side wind direction (correlated)
## Wind Speed
#WS_sPCAN = np.array(PCAN['Starboard True Wind Speed (knot)']) # starboard side wind speed (correlated)
#WS_sPCAN = WS_sPCAN * 0.514444444 # Convert from knots to m/s
#WS_pPCAN = np.array(PCAN['Port True Wind Speed (knot)']) # port side wind speed (correlated)
#WS_pPCAN = WS_pPCAN * 0.514444444 # Convert from knots to m/s
#WS_PCAN = (WS_sPCAN + WS_pPCAN)/2 #Average the wind speed for port and starboard
## Vector Mean Wind Direction
#WD_vectPCAN = ((WD_sPCAN * WS_sPCAN) / (WS_sPCAN + WS_pPCAN)) + ((WD_pPCAN * WS_pPCAN) / (WS_sPCAN + WS_pPCAN)) # Calculate the vector mean wind direction
#
# SIPEXII (2012)----
# Wind Direction
WD_sSIP = np.array(SIPEXII['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_pSIP = np.array(SIPEXII['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
#WD_sSIP = np.array(SIPEXII['wnd_dir_strbrd_uncorr_deg']) # starboard side wind direction (uncorrelated)
#WD_pSIP = np.array(SIPEXII['wnd_dir_port_uncorr_deg']) # port side wind direction (uncorrelated)
# Average wind direction
WD_Avg_SIP = (WD_sSIP + WD_pSIP) / 2

# Wind Speed
WS_sSIP = np.array(SIPEXII['wnd_spd_strbrd_corr_knot']) * 0.514444444 # starboard side wind speed (correlated)
WS_pSIP = np.array(SIPEXII['wnd_spd_port_corr_knot'])   * 0.514444444 # port side wind speed (correlated)
#WS_sSIP = np.array(SIPEXII['wnd_spd_strbrd_uncorr_knot']) * 0.514444444 # starboard side wind speed (uncorrelated)
#WS_pSIP = np.array(SIPEXII['wnd_spd_port_uncorr_knot'])   * 0.514444444 # port side wind speed (uncorrelated)
# Average wind speed
WS_Avg_SIP = (WS_sSIP + WS_pSIP) / 2

# Vector Mean Wind Direction
WD_vectSIP = ((WD_sSIP * WS_sSIP) / (WS_sSIP + WS_pSIP)) + ((WD_pSIP * WS_pSIP) / (WS_sSIP + WS_pSIP)) # Calculate the vector mean wind direction

# ONLY ABOVE 60S

## PCAN (2017)--------
## Wind Direction
#WD_sPCAN_N60S = np.array(PCAN_N60S['Starboard True Wind Direction (degree)']) # starboard side wind direction (correlated)
#WD_pPCAN_N60S = np.array(PCAN_N60S['Port True Wind Direction (degree)']) # port side wind direction (correlated)
## Wind Speed
#WS_sPCAN_N60S = np.array(PCAN_N60S['Starboard True Wind Speed (knot)']) # starboard side wind speed (correlated)
#WS_sPCAN_N60S = WS_sPCAN_N60S * 0.514444444 # Convert from knots to m/s
#WS_pPCAN_N60S = np.array(PCAN_N60S['Port True Wind Speed (knot)']) # port side wind speed (correlated)
#WS_pPCAN_N60S = WS_pPCAN_N60S * 0.514444444 # Convert from knots to m/s
#WS_PCAN_N60S = (WS_sPCAN_N60S + WS_pPCAN_N60S)/2 #Average the wind speed for port and starboard
## Vector Mean Wind Direction
#WD_vectPCAN_N60S = ((WD_sPCAN_N60S * WS_sPCAN_N60S) / (WS_sPCAN_N60S + WS_pPCAN_N60S)) + ((WD_pPCAN_N60S * WS_pPCAN_N60S) / (WS_sPCAN_N60S + WS_pPCAN_N60S)) # Calculate the vector mean wind direction

# SIPEXII (2012)----
# Wind Direction
WD_sSIP_N60S = np.array(SIPEXII_N60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_pSIP_N60S = np.array(SIPEXII_N60S['wnd_dir_port_corr_deg'])   # port side wind direction (correlated)
#WD_sSIP_N60S = np.array(SIPEXII_N60S['wnd_dir_strbrd_uncorr_deg']) # starboard side wind direction (uncorrelated)
#WD_pSIP_N60S = np.array(SIPEXII_N60S['wnd_dir_port_uncorr_deg'])   # port side wind direction (uncorrelated)
# Average wind direction
WD_Avg_SIP_N60S = (WD_sSIP_N60S + WD_pSIP_N60S) / 2

# Wind Speed
WS_sSIP_N60S = np.array(SIPEXII_N60S['wnd_spd_strbrd_corr_knot']) * 0.514444444 # starboard side wind speed (correlated)
WS_pSIP_N60S = np.array(SIPEXII_N60S['wnd_spd_port_corr_knot'])   * 0.514444444 # port side wind speed (correlated)
#WS_sSIP_N60S = np.array(SIPEXII_N60S['wnd_spd_strbrd_uncorr_knot']) * 0.514444444 # starboard side wind speed (uncorrelated)
#WS_pSIP_N60S = np.array(SIPEXII_N60S['wnd_spd_port_uncorr_knot'])   * 0.514444444 # port side wind speed (uncorrelated)
# Average wind speed
WS_Avg_SIP_N60S = (WS_sSIP_N60S + WS_pSIP_N60S) / 2

# Vector Mean Wind Direction
WD_vectSIP_N60S = ((WD_sSIP_N60S * WS_sSIP_N60S) / (WS_sSIP_N60S + WS_pSIP_N60S)) + ((WD_pSIP_N60S * WS_pSIP_N60S) / (WS_sSIP_N60S + WS_pSIP_N60S)) # Calculate the vector mean wind direction

# ONLY BELOW 60S

## PCAN (2017)--------
## Wind Direction
#WD_sPCAN_60S = np.array(PCAN_60S['Starboard True Wind Direction (degree)']) # starboard side wind direction (correlated)
#WD_pPCAN_60S = np.array(PCAN_60S['Port True Wind Direction (degree)']) # port side wind direction (correlated)
## Wind Speed
#WS_sPCAN_60S = np.array(PCAN_60S['Starboard True Wind Speed (knot)']) # starboard side wind speed (correlated)
#WS_sPCAN_60S = WS_sPCAN_60S * 0.514444444 # Convert from knots to m/s
#WS_pPCAN_60S = np.array(PCAN_60S['Port True Wind Speed (knot)']) # port side wind speed (correlated)
#WS_pPCAN_60S = WS_pPCAN_60S * 0.514444444 # Convert from knots to m/s
#WS_PCAN_60S = (WS_sPCAN_60S + WS_pPCAN_60S)/2 #Average the wind speed for port and starboard
## Vector Mean Wind Direction
#WD_vectPCAN_60S = ((WD_sPCAN_60S * WS_sPCAN_60S) / (WS_sPCAN_60S + WS_pPCAN_60S)) + ((WD_pPCAN_60S * WS_pPCAN_60S) / (WS_sPCAN_60S + WS_pPCAN_60S)) # Calculate the vector mean wind direction

# SIPEXII (2012)----
# Wind Direction
WD_sSIP_60S = np.array(SIPEXII_60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_pSIP_60S = np.array(SIPEXII_60S['wnd_dir_port_corr_deg'])   # port side wind direction (correlated)
#WD_sSIP_60S = np.array(SIPEXII_60S['wnd_dir_strbrd_uncorr_deg']) # starboard side wind direction (uncorrelated)
#WD_pSIP_60S = np.array(SIPEXII_60S['wnd_dir_port_uncorr_deg'])   # port side wind direction (uncorrelated)
# Average wind direction
WD_Avg_SIP_60S = (WD_sSIP_60S + WD_pSIP_60S) / 2

# Wind Speed
WS_sSIP_60S = np.array(SIPEXII_60S['wnd_spd_strbrd_corr_knot']) * 0.514444444 # starboard side wind speed (correlated)
WS_pSIP_60S = np.array(SIPEXII_60S['wnd_spd_port_corr_knot'])   * 0.514444444 # port side wind speed (correlated)
#WS_sSIP_60S = np.array(SIPEXII_60S['wnd_spd_strbrd_uncorr_knot']) * 0.514444444 # starboard side wind speed (uncorrelated)
#WS_pSIP_60S = np.array(SIPEXII_60S['wnd_spd_port_uncorr_knot'])   * 0.514444444 # port side wind speed (uncorrelated)
# Average wind speed
WS_Avg_SIP_60S = (WS_sSIP_60S + WS_pSIP_60S) / 2

# Vector Mean Wind Direction
WD_vectSIP_60S = ((WD_sSIP_60S * WS_sSIP_60S) / (WS_sSIP_60S + WS_pSIP_60S)) + ((WD_pSIP_60S * WS_pSIP_60S) / (WS_sSIP_60S + WS_pSIP_60S)) # Calculate the vector mean wind direction

#------------------------------------------------------------------------------
# PLOT A WIND ROSE OF WIND SPEED AND VECTOR MEAN WIND DIRECTION
fig = plt.figure(figsize=(20,15)) # (width, height)

#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vectPCAN, WS_PCAN, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("PCAN (all latitudes)", position=(0.5, 1.1))

#rect_2 = [0.175, 0.7, 0.225, 0.225]     # [left, bottom, width, height]
#ax = WindroseAxes(fig, rect_2)
#fig.add_axes(ax)
#ax.bar(WD_vectSIP, WS_SIP, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("SIPEXII (all latitudes)", position=(0.5, 1.1))

#rect_3 = [0, 0.375, 0.225, 0.225]     # [left, bottom, width, height]
#ax = WindroseAxes(fig, rect_3)
#fig.add_axes(ax)
#ax.bar(WD_vectPCAN_N60S, WS_PCAN_N60S, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("PCAN (above 60S)", position=(0.5, 1.1))

rect_4 = [0.175, 0.375, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
ax = WindroseAxes(fig, rect_4)    # creates the axes of specified dimensions 
fig.add_axes(ax)
ax.bar(WD_Avg_SIP_N60S, WS_Avg_SIP_N60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
ax.set_title("SIPEXII", position=(0.5, 1.1), fontsize=20, pad=20)
ax.set_ylabel('North of 64$^\circ$S', fontsize=20, labelpad=30)

#rect_5 = [0, 0.05, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_5)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vectPCAN_60S, WS_PCAN_60S, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("PCAN (below 60S)", position=(0.5, 1.1))

rect_6 = [0.175, 0.05, 0.225, 0.225]     # [left, bottom, width, height]
ax = WindroseAxes(fig, rect_6)
fig.add_axes(ax)
ax.bar(WD_Avg_SIP_60S, WS_Avg_SIP_60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
#ax.set_title("SIPEXII", position=(0.5, 1.1), fontsize=20, pad=20)
ax.set_ylabel('South of 65$^\circ$S', fontsize=20, labelpad=30)

ax.set_legend()
ax.legend(title="wind speed (m/s)", loc=(1.2, 0))
