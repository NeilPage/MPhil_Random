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
V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/ShipTrack/V1_18_underway_60.csv', index_col=0)
V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/ShipTrack/V2_18_underway_60.csv', index_col=0)
V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/ShipTrack/V3_18_underway_60.csv', index_col=0)
#V4_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V04/CAMMPCAN_V4_underway_60.csv') # V4 Cruise data

# fig = plt.figure()
# ax = plt.subplot(111)
# ax2 = ax.twinx()
# ax.plot(V1_17.index, V1_17['latitude'],      marker='o', c='r', markersize = 2.0, linestyle='none')
# ax2.plot(V1_17.index, V1_17['atm_press_hpa'], marker='o', c='b', markersize = 2.0, linestyle='none')

#------------------------------------------------------------------------------
# FILTER OUT DATA CORRESPONDING TO BELOW 65 DEGREES SOUTH

# V1
LATS = V1_17['latitude']<=-65
V1_17_60S = V1_17[LATS]
# V2
LATS = V2_17['latitude']<=-65
V2_17_60S = V2_17[LATS]
# V3
LATS = V3_17['latitude']<=-65
V3_17_60S = V3_17[LATS]

# FILTER OUT DATA CORRESPONDING TO ABOVE 65 DEGREES SOUTH

# V1
LATS = V1_17['latitude']>-65
V1_17_N60S = V1_17[LATS]
# V2
LATS = V2_17['latitude']>-65
V2_17_N60S = V2_17[LATS]
# V3
LATS = V3_17['latitude']>-65
V3_17_N60S = V3_17[LATS]

#------------------------------------------------------------------------------
# Define the variables

# ALL LATITUDES
# V1----------------
# Wind Direction
WD_s1_17 = np.array(V1_17['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p1_17 = np.array(V1_17['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V1_17 = (WD_s1_17 + WD_p1_17) / 2

# Wind Speed
WS_s1_17 = np.array(V1_17['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s1_17 = WS_s1_17 * 0.514444444 # Convert from knots to m/s
WS_p1_17 = np.array(V1_17['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p1_17 = WS_p1_17 * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V1_17 = (WS_s1_17 + WS_p1_17) / 2

# Vector Mean Wind Direction
WD_vect1_17 = ((WD_s1_17 * WS_s1_17) / (WS_s1_17 + WS_p1_17)) + ((WD_p1_17 * WS_p1_17) / (WS_s1_17 + WS_p1_17)) # Calculate the vector mean wind direction

# V2----------------
# Wind Direction
WD_s2_17 = np.array(V2_17['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p2_17 = np.array(V2_17['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V2_17 = (WD_s2_17 + WD_p2_17) / 2

# Wind Speed
WS_s2_17 = np.array(V2_17['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s2_17 = WS_s2_17 * 0.514444444 # Convert from knots to m/s
WS_p2_17 = np.array(V2_17['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p2_17 = WS_p2_17 * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V2_17 = (WS_s2_17 + WS_p2_17) / 2

# Vector Mean Wind Direction
WD_vect2_17 = ((WD_s2_17 * WS_s2_17) / (WS_s2_17 + WS_p2_17)) + ((WD_p2_17 * WS_p2_17) / (WS_s2_17 + WS_p2_17)) # Calculate the vector mean wind direction

# V3----------------
# Wind Direction
WD_s3_17 = np.array(V3_17['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p3_17 = np.array(V3_17['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V3_17 = (WD_s3_17 + WD_p3_17) / 2

# Wind Speed
WS_s3_17 = np.array(V3_17['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s3_17 = WS_s3_17 * 0.514444444 # Convert from knots to m/s
WS_p3_17 = np.array(V3_17['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p3_17 = WS_p3_17 * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V3_17 = (WS_s3_17 + WS_p3_17) / 2

# Vector Mean Wind Direction
WD_vect3_17 = ((WD_s3_17 * WS_s3_17) / (WS_s3_17 + WS_p3_17)) + ((WD_p3_17 * WS_p3_17) / (WS_s3_17 + WS_p3_17)) # Calculate the vector mean wind direction

## V4----------------
## Wind Direction
#WD_s4_17 = np.array(V4_17['WND_DIR_STRBD_CORR_DEG']) # starboard side wind direction (correlated)
#WD_p4_17 = np.array(V4_17['WND_DIR_PORT_CORR_DEG']) # port side wind direction (correlated)
## Wind Speed
#WS_s4_17 = np.array(V4_17['WND_SPD_STRBD_CORR_KNOT']) # starboard side wind speed (correlated)
#WS_s4_17 = WS_s4_17 * 0.514444444 # Convert from knots to m/s
#WS_p4_17 = np.array(V4_17['WND_SPD_PORT_CORR_KNOT']) # port side wind speed (correlated)
#WS_p4_17 = WS_p4_17 * 0.514444444 # Convert from knots to m/s
#WS_4_17 = (WS_s4_17 + WS_p4_17)/2 #Average the wind speed for port and starboard
## Vector Mean Wind Direction
#WD_vect4_17 = ((WD_s4_17 * WS_s4_17) / (WS_s4_17 + WS_p4_17)) + ((WD_p4_17 * WS_p4_17) / (WS_s4_17 + WS_p4_17)) # Calculate the vector mean wind direction

# ONLY ABOVE 60S
# V1----------------
# Wind Direction
WD_s1_17_N60S = np.array(V1_17_N60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p1_17_N60S = np.array(V1_17_N60S['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V1_17_N60S = (WD_s1_17_N60S + WD_p1_17_N60S) / 2

# Wind Speed
WS_s1_17_N60S = np.array(V1_17_N60S['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s1_17_N60S = WS_s1_17_N60S * 0.514444444 # Convert from knots to m/s
WS_p1_17_N60S = np.array(V1_17_N60S['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p1_17_N60S = WS_p1_17_N60S * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V1_17_N60S = (WS_s1_17_N60S + WS_p1_17_N60S) / 2

# Vector Mean Wind Direction
WD_vect1_17_N60S = ((WD_s1_17_N60S * WS_s1_17_N60S) / (WS_s1_17_N60S + WS_p1_17_N60S)) + ((WD_p1_17_N60S * WS_p1_17_N60S) / (WS_s1_17_N60S + WS_p1_17_N60S)) # Calculate the vector mean wind direction

# V2----------------
# Wind Direction
WD_s2_17_N60S = np.array(V2_17_N60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p2_17_N60S = np.array(V2_17_N60S['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V2_17_N60S = (WD_s2_17_N60S + WD_p2_17_N60S) / 2

# Wind Speed
WS_s2_17_N60S = np.array(V2_17_N60S['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s2_17_N60S = WS_s2_17_N60S * 0.514444444 # Convert from knots to m/s
WS_p2_17_N60S = np.array(V2_17_N60S['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p2_17_N60S = WS_p2_17_N60S * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V2_17_N60S = (WS_s2_17_N60S + WS_p2_17_N60S) / 2

# Vector Mean Wind Direction
WD_vect2_17_N60S = ((WD_s2_17_N60S * WS_s2_17_N60S) / (WS_s2_17_N60S + WS_p2_17_N60S)) + ((WD_p2_17_N60S * WS_p2_17_N60S) / (WS_s2_17_N60S + WS_p2_17_N60S)) # Calculate the vector mean wind direction

# V3----------------
# Wind Direction
WD_s3_17_N60S = np.array(V3_17_N60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p3_17_N60S = np.array(V3_17_N60S['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V3_17_N60S = (WD_s3_17_N60S + WD_p3_17_N60S) / 2

# Wind Speed
WS_s3_17_N60S = np.array(V3_17_N60S['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s3_17_N60S = WS_s3_17_N60S * 0.514444444 # Convert from knots to m/s
WS_p3_17_N60S = np.array(V3_17_N60S['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p3_17_N60S = WS_p3_17_N60S * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V3_17_N60S = (WS_s3_17_N60S + WS_p3_17_N60S) / 2

# Vector Mean Wind Direction
WD_vect3_17_N60S = ((WD_s3_17_N60S * WS_s3_17_N60S) / (WS_s3_17_N60S + WS_p3_17_N60S)) + ((WD_p3_17_N60S * WS_p3_17_N60S) / (WS_s3_17_N60S + WS_p3_17_N60S)) # Calculate the vector mean wind direction


# ONLY BELOW 60S
# V1----------------
# Wind Direction
WD_s1_17_60S = np.array(V1_17_60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p1_17_60S = np.array(V1_17_60S['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V1_17_60S = (WD_s1_17_60S + WD_p1_17_60S) / 2

# Wind Speed
WS_s1_17_60S = np.array(V1_17_60S['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s1_17_60S = WS_s1_17_60S * 0.514444444 # Convert from knots to m/s
WS_p1_17_60S = np.array(V1_17_60S['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p1_17_60S = WS_p1_17_60S * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V1_17_60S = (WS_s1_17_60S + WS_p1_17_60S) / 2

# Vector Mean Wind Direction
WD_vect1_17_60S = ((WD_s1_17_60S * WS_s1_17_60S) / (WS_s1_17_60S + WS_p1_17_60S)) + ((WD_p1_17_60S * WS_p1_17_60S) / (WS_s1_17_60S + WS_p1_17_60S)) # Calculate the vector mean wind direction

# V2----------------
# Wind Direction
WD_s2_17_60S = np.array(V2_17_60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p2_17_60S = np.array(V2_17_60S['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V2_17_60S = (WD_s2_17_60S + WD_p2_17_60S) / 2

# Wind Speed
WS_s2_17_60S = np.array(V2_17_60S['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s2_17_60S = WS_s2_17_60S * 0.514444444 # Convert from knots to m/s
WS_p2_17_60S = np.array(V2_17_60S['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p2_17_60S = WS_p2_17_60S * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V2_17_60S = (WS_s2_17_60S + WS_p2_17_60S) / 2

# Vector Mean Wind Direction
WD_vect2_17_60S = ((WD_s2_17_60S * WS_s2_17_60S) / (WS_s2_17_60S + WS_p2_17_60S)) + ((WD_p2_17_60S * WS_p2_17_60S) / (WS_s2_17_60S + WS_p2_17_60S)) # Calculate the vector mean wind direction

# V3----------------
# Wind Direction
WD_s3_17_60S = np.array(V3_17_60S['wnd_dir_strbrd_corr_deg']) # starboard side wind direction (correlated)
WD_p3_17_60S = np.array(V3_17_60S['wnd_dir_port_corr_deg']) # port side wind direction (correlated)
# Average wind direction
WD_Avg_V3_17_60S = (WD_s3_17_60S + WD_p3_17_60S) / 2

# Wind Speed
WS_s3_17_60S = np.array(V3_17_60S['wnd_spd_strbrd_corr_knot']) # starboard side wind speed (correlated)
WS_s3_17_60S = WS_s3_17_60S * 0.514444444 # Convert from knots to m/s
WS_p3_17_60S = np.array(V3_17_60S['wnd_spd_port_corr_knot']) # port side wind speed (correlated)
WS_p3_17_60S = WS_p3_17_60S * 0.514444444 # Convert from knots to m/s
# Average wind speed
WS_Avg_V3_17_60S = (WS_s3_17_60S + WS_p3_17_60S) / 2

# Vector Mean Wind Direction
WD_vect3_17_60S = ((WD_s3_17_60S * WS_s3_17_60S) / (WS_s3_17_60S + WS_p3_17_60S)) + ((WD_p3_17_60S * WS_p3_17_60S) / (WS_s3_17_60S + WS_p3_17_60S)) # Calculate the vector mean wind direction

#------------------------------------------------------------------------------
# PLOT A WIND ROSE OF WIND SPEED AND VECTOR MEAN WIND DIRECTION
fig = plt.figure(figsize=(20,15)) # (width, height)

#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect1_17, WS_1_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("V1 (all latitudes)", position=(0.5, 1.1))
#
#rect_2 = [0.175, 0.7, 0.225, 0.225]     # [left, bottom, width, height]
#ax = WindroseAxes(fig, rect_2)
#fig.add_axes(ax)
#ax.bar(WD_vect2_17, WS_2_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("V2 (all latitudes)", position=(0.5, 1.1))
#
#rect_3 = [0.35, 0.7, 0.225, 0.225]     # [left, bottom, width, height]
#ax = WindroseAxes(fig, rect_3)
#fig.add_axes(ax)
#ax.bar(WD_vect3_17, WS_3_17, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,35, 5))
#ax.set_title("V3 (all latitudes)", position=(0.5, 1.1))

rect_4 = [0, 0.375, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
ax = WindroseAxes(fig, rect_4)    # creates the axes of specified dimensions 
fig.add_axes(ax)
ax.bar(WD_Avg_V1_17_N60S, WS_Avg_V1_17_N60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
ax.set_title("V1", position=(0.5, 1.1), fontsize=20, pad=20)
ax.set_ylabel('North of 65$^\circ$S', fontsize=20, labelpad=30)

rect_5 = [0.175, 0.375, 0.225, 0.225]     # [left, bottom, width, height]
ax = WindroseAxes(fig, rect_5)
fig.add_axes(ax)
ax.bar(WD_Avg_V2_17_N60S, WS_Avg_V2_17_N60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
ax.set_title("V2", position=(0.5, 1.1), fontsize=20, pad=20)

rect_6= [0.35, 0.375, 0.225, 0.225]     # [left, bottom, width, height]
ax = WindroseAxes(fig, rect_6)
fig.add_axes(ax)
ax.bar(WD_Avg_V3_17_N60S, WS_Avg_V3_17_N60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
ax.set_title("V3", position=(0.5, 1.1), fontsize=20, pad=20)

rect_7 = [0, 0.05, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
ax = WindroseAxes(fig, rect_7)    # creates the axes of specified dimensions 
fig.add_axes(ax)
ax.bar(WD_Avg_V1_17_60S, WS_Avg_V1_17_60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
#ax.set_title("V1", position=(0.5, 1.1), fontsize=20)
ax.set_ylabel('South of 65$^\circ$S', fontsize=20, labelpad=30)

rect_8 = [0.175, 0.05, 0.225, 0.225]     # [left, bottom, width, height]
ax = WindroseAxes(fig, rect_8)
fig.add_axes(ax)
ax.bar(WD_Avg_V2_17_60S, WS_Avg_V2_17_60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
#ax.set_title("V2", position=(0.5, 1.1), fontsize=20)

rect_9 = [0.35, 0.05, 0.225, 0.225]     # [left, bottom, width, height]
ax = WindroseAxes(fig, rect_9)
fig.add_axes(ax)
ax.bar(WD_Avg_V3_17_60S, WS_Avg_V3_17_60S, normed=True, opening=1.0, edgecolor='black', bins=np.arange(0, 25.0, 2.5))
#ax.set_title("V3", position=(0.5, 1.1), fontsize=20)

ax.set_legend()
ax.legend(title="wind speed (m/s)", loc=(1.2, 0))
