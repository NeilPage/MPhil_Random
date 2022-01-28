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

Event1 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V02/Event1.csv') # V1 Cruise data
Event2 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V03/Event2.csv') # V2 Cruise data 
Event3 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V03/Event3.csv') # V3 Cruise data
Event4 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V01/Event4.csv') # V1 Cruise data
Event5 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V03/Event5.csv') # V2 Cruise data 
Event6 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V03/Event6.csv') # V3 Cruise data
Event7 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/Event7.csv') # SIPEXII Cruise data

#------------------------------------------------------------------------------
# Define the variables


# EVENT1----------------
# Wind Direction
WD_s1 = np.array(Event1['WND_DIR_STRBD_CORR_DEG'])
WD_p1 = np.array(Event1['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s1 = np.array(Event1['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p1 = np.array(Event1['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_1 = (WS_s1 + WS_p1)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect1 = ((WD_s1 * WS_s1) / (WS_s1 + WS_p1)) + ((WD_p1 * WS_p1) / (WS_s1 + WS_p1)) # Calculate the vector mean wind direction

# EVENT2----------------
# Wind Direction
WD_s2 = np.array(Event2['WND_DIR_STRBD_CORR_DEG'])
WD_p2 = np.array(Event2['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s2 = np.array(Event2['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p2 = np.array(Event2['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_2 = (WS_s2 + WS_p2)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect2 = ((WD_s2 * WS_s2) / (WS_s2 + WS_p2)) + ((WD_p2 * WS_p2) / (WS_s2 + WS_p2)) # Calculate the vector mean wind direction

# EVENT3----------------
# Wind Direction
WD_s3 = np.array(Event3['WND_DIR_STRBD_CORR_DEG'])
WD_p3 = np.array(Event3['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s3 = np.array(Event3['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p3 = np.array(Event3['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_3 = (WS_s3 + WS_p3)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect3 = ((WD_s3 * WS_s3) / (WS_s3 + WS_p3)) + ((WD_p3 * WS_p3) / (WS_s3 + WS_p3)) # Calculate the vector mean wind direction

# EVENT4----------------
# Wind Direction
WD_s4 = np.array(Event4['WND_DIR_STRBD_CORR_DEG'])
WD_p4 = np.array(Event4['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s4 = np.array(Event4['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p4 = np.array(Event4['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_4 = (WS_s4 + WS_p4)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect4 = ((WD_s4 * WS_s4) / (WS_s4 + WS_p4)) + ((WD_p4 * WS_p4) / (WS_s4 + WS_p4)) # Calculate the vector mean wind direction

# EVENT5----------------
# Wind Direction
WD_s5 = np.array(Event5['WND_DIR_STRBD_CORR_DEG'])
WD_p5 = np.array(Event5['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s5 = np.array(Event5['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p5 = np.array(Event5['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_5 = (WS_s5 + WS_p5)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect5 = ((WD_s5 * WS_s5) / (WS_s5 + WS_p5)) + ((WD_p5 * WS_p5) / (WS_s5 + WS_p5)) # Calculate the vector mean wind direction

# EVENT6----------------
# Wind Direction
WD_s6 = np.array(Event6['WND_DIR_STRBD_CORR_DEG'])
WD_p6 = np.array(Event6['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s6 = np.array(Event6['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p6 = np.array(Event6['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_6 = (WS_s6 + WS_p6)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect6 = ((WD_s6 * WS_s6) / (WS_s6 + WS_p6)) + ((WD_p6 * WS_p6) / (WS_s6 + WS_p6)) # Calculate the vector mean wind direction

# EVENT7----------------
# Wind Direction
WD_s7 = np.array(Event7['WND_DIR_STRBD_CORR_DEG'])
WD_p7 = np.array(Event7['WND_DIR_PORT_CORR_DEG'])
# Wind Speed
WS_s7 = np.array(Event7['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_p7 = np.array(Event7['WND_SPD_PORT_CORR_KNOT']) * 0.514444444 # Convert from knots to m/s
WS_7 = (WS_s7 + WS_p7)/2 #Average the wind speed for port and starboard
# Vector Mean Wind Direction
WD_vect7 = ((WD_s7 * WS_s7) / (WS_s7 + WS_p7)) + ((WD_p7 * WS_p7) / (WS_s7 + WS_p7)) # Calculate the vector mean wind direction

#------------------------------------------------------------------------------
# PLOT A WIND ROSE OF WIND SPEED AND VECTOR MEAN WIND DIRECTION
fig = plt.figure(figsize=(20,15)) # (width, height)

# Event1
rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
fig.add_axes(ax)
ax.bar(WD_vect1, WS_1, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
ax.set_title("Wind speed and direction", position=(0.5, 1.1))

## Event2
#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect2, WS_2, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
#ax.set_title("Wind speed and direction", position=(0.5, 1.1))
#
##Event3
#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect3, WS_3, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
#ax.set_title("Wind speed and direction", position=(0.5, 1.1))
#
##Event4
#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect4, WS_4, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
#ax.set_title("Wind speed and direction", position=(0.5, 1.1))
#
##Event5
#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect5, WS_5, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
#ax.set_title("Wind speed and direction", position=(0.5, 1.1))

##Event6
#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect6, WS_6, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
#ax.set_title("Wind speed and direction", position=(0.5, 1.1))

##Event7
#rect_1 = [0, 0.7, 0.225, 0.225]     # [left, bottom, width, height] as a fraction of total figure size
#ax = WindroseAxes(fig, rect_1)    # creates the axes of specified dimensions 
#fig.add_axes(ax)
#ax.bar(WD_vect7, WS_7, normed=True, opening=0.8, edgecolor='white', bins=np.arange(0,25, 3))
#ax.set_title("Wind speed and direction", position=(0.5, 1.1))

ax.set_legend()
ax.legend(title="wind speed (m/s)", loc=(1.2, 0))
