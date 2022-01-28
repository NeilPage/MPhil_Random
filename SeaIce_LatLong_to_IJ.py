#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:23:56 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
import xarray as xr

# DATA HANDLING PACKAGES
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# DRAWING PACKAGES
import cartopy.crs as ccrs
#from matplotlib import cm                   # imports the colormap function from matplotlib
import matplotlib.ticker as mticker
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable

#------------------------------------------------------------------------------
# DEFINE THE DATASET

## CAMMPCAN 2017-18
#V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V1_17_Minute_LatLong.csv')
#V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V2_17_Minute_LatLong.csv')
#V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V3_17_Minute_LatLong.csv')
#
## CAMMPCAN 2018-19
#V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1_18_Minute_LatLong.csv')
#V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V2_18_Minute_LatLong.csv')
#V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V3_18_Minute_LatLong.csv')

# SIPEXII 2012
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/SIPEXII_Minute_LatLong.csv')

#------------------------------------------------------------------------------
# DEFINE THE NSIDC VARIABLES

## Date, Lat and Long
#Date_V1_17 = np.array(V1_17['Date_V1_17'])
#Lat_V1_17 = np.array(V1_17['Lat_V1_17'])
#Lon_V1_17 = np.array(V1_17['Lon_V1_17'])
#
#Date_V2_17 = np.array(V2_17['Date_V2_17'])
#Lat_V2_17 = np.array(V2_17['Lat_V2_17'])
#Lon_V2_17 = np.array(V2_17['Lon_V2_17'])
#
#Date_V3_17 = np.array(V3_17['Date_V3_17'])
#Lat_V3_17 = np.array(V3_17['Lat_V3_17'])
#Lon_V3_17 = np.array(V3_17['Lon_V3_17'])
#
#Date_V1_18 = np.array(V1_18['Date_V1_18'])
#Lat_V1_18 = np.array(V1_18['Lat_V1_18'])
#Lon_V1_18 = np.array(V1_18['Lon_V1_18'])
#
#Date_V2_18 = np.array(V2_18['Date_V2_18'])
#Lat_V2_18 = np.array(V2_18['Lat_V2_18'])
#Lon_V2_18 = np.array(V2_18['Lon_V2_18'])
#
#Date_V3_18 = np.array(V3_18['Date_V3_18'])
#Lat_V3_18 = np.array(V3_18['Lat_V3_18'])
#Lon_V3_18 = np.array(V3_18['Lon_V3_18'])

Date_SIPEXII = np.array(SIPEXII['Date_SIPEXII'])
Lat_SIPEXII = np.array(SIPEXII['Lat_SIPEXII'])
Lon_SIPEXII = np.array(SIPEXII['Lon_SIPEXII'])

# Hemisphere and grid cell size
hemisphere = -1 # 1 for Northern hemisphere, -1 for Southern
grid = 25 # 6.25, 12.5 or 25; the grid cell dimensions in km

#------------------------------------------------------------------------------
# FUNCTION TO TRANSFORM FROM LAT/LONG COORDINATES TO NSIDC POLAR STEREOGRAPHIC (I,J)

def polar_lonlat_to_xy(longitude, latitude, true_scale_lat, re, e, hemisphere):
    """Convert from geodetic longitude and latitude to Polar Stereographic
    (X, Y) coordinates in km.
    Args:
        longitude (float): longitude or longitude array in degrees
        latitude (float): latitude or latitude array in degrees (positive)
        true_scale_lat (float): true-scale latitude in degrees
        re (float): Earth radius in km
        e (float): Earth eccentricity
        hemisphere (1 or -1): Northern or Southern hemisphere
    Returns:
        If longitude and latitude are scalars then the result is a
        two-element list containing [X, Y] in km.
        If longitude and latitude are numpy arrays then the result will be a
        two-element list where the first element is a numpy array containing
        the X coordinates and the second element is a numpy array containing
        the Y coordinates.
    """

    lat = abs(latitude) * np.pi / 180
    lon = longitude * np.pi / 180
    slat = true_scale_lat * np.pi / 180

    e2 = e * e

    # Snyder (1987) p. 161 Eqn 15-9
    t = np.tan(np.pi / 4 - lat / 2) / \
        ((1 - e * np.sin(lat)) / (1 + e * np.sin(lat))) ** (e / 2)

    if abs(90 - true_scale_lat) < 1e-5:
        # Snyder (1987) p. 161 Eqn 21-33
        rho = 2 * re * t / np.sqrt((1 + e) ** (1 + e) * (1 - e) ** (1 - e))
    else:
        # Snyder (1987) p. 161 Eqn 21-34
        tc = np.tan(np.pi / 4 - slat / 2) / \
            ((1 - e * np.sin(slat)) / (1 + e * np.sin(slat))) ** (e / 2)
        mc = np.cos(slat) / np.sqrt(1 - e2 * (np.sin(slat) ** 2))
        rho = re * mc * t / tc

    x = rho * hemisphere * np.sin(hemisphere * lon)
    y = -rho * hemisphere * np.cos(hemisphere * lon)
    return [x, y]  

def nsidc_polar_lonlat(longitude, latitude, grid, hemisphere):
    """Transform from geodetic longitude and latitude coordinates
    to NSIDC Polar Stereographic I, J coordinates
    
    Args:
        longitude (float): longitude or longitude array in degrees
        latitude (float): latitude or latitude array in degrees (positive)
        grid (float): 6.25, 12.5 or 25; the grid cell dimensions in km
        hemisphere (1 or -1): Northern or Southern hemisphere
    
    Returns:
        If longitude and latitude are scalars then the result is a
        two-element list containing [I, J].
        If longitude and latitude are numpy arrays then the result will
        be a two-element list where the first element is a numpy array for
        the I coordinates and the second element is a numpy array for
        the J coordinates.
    Examples:
        print(nsidc_polar_lonlat(350.0, 34.41, 12.5, 1))
            [608, 896]
    """

    true_scale_lat = 70
    re = 6378.273
    e = 0.081816153

    if grid != 6.25 and grid != 12.5 and grid != 25:
        raise ValueError("Legal grid value are 6.25, 12.5, or 25")
    
    if hemisphere >= 0:
        delta = 45
        imax = 1216
        jmax = 1792
        xmin = -3850 + grid/2
        ymin = -5350 + grid/2
    else:
        delta = 0
        imax = 1264
        jmax = 1328
        xmin = -3950 + grid/2
        ymin = -3950 + grid/2

    if grid == 12.5:
        imax = imax//2
        jmax = jmax//2
    elif grid == 25:
        imax = imax//4
        jmax = jmax//4

    xy = polar_lonlat_to_xy(longitude + delta, np.abs(latitude),
                            true_scale_lat, re, e, hemisphere)
    i = (np.round((xy[0] - xmin)/grid)).astype(int) + 1
    j = (np.round((xy[1] - ymin)/grid)).astype(int) + 1
    # Flip grid orientation in the 'y' direction
    j = jmax - j + 1
    return [i, j]

#------------------------------------------------------------------------------
# CALL THE FUNCTION

## Set up reservoir arrays
#iV1_17 = np.zeros(len(Lat_V1_17))
#jV1_17 = np.zeros(len(Lat_V1_17))
#iV2_17 = np.zeros(len(Lat_V2_17))
#jV2_17 = np.zeros(len(Lat_V2_17))
#iV3_17 = np.zeros(len(Lat_V3_17))
#jV3_17 = np.zeros(len(Lat_V3_17))
#iV1_18 = np.zeros(len(Lat_V1_18))
#jV1_18 = np.zeros(len(Lat_V1_18))
#iV2_18 = np.zeros(len(Lat_V2_18))
#jV2_18 = np.zeros(len(Lat_V2_18))
#iV3_18 = np.zeros(len(Lat_V3_18))
#jV3_18 = np.zeros(len(Lat_V3_18))
iSIPEXII = np.zeros(len(Lat_SIPEXII))
jSIPEXII = np.zeros(len(Lat_SIPEXII))

## Use the function to loop over each array
#for i in range(len(Lat_V1_17)):
#    iV1_17[i],jV1_17[i] = nsidc_polar_lonlat(Lon_V1_17[i],Lat_V1_17[i],grid,hemisphere) # V1_17
#
#for i in range(len(Lat_V2_17)):
#    iV2_17[i],jV2_17[i] = nsidc_polar_lonlat(Lon_V2_17[i],Lat_V2_17[i],grid,hemisphere) # V2_17
#
#for i in range(len(Lat_V3_17)):
#    iV3_17[i],jV3_17[i] = nsidc_polar_lonlat(Lon_V3_17[i],Lat_V3_17[i],grid,hemisphere) # V1_17
#
#for i in range(len(Lat_V1_18)):
#    iV1_18[i],jV1_18[i] = nsidc_polar_lonlat(Lon_V1_18[i],Lat_V1_18[i],grid,hemisphere) # V1_17
#
#for i in range(len(Lat_V2_18)):
#    iV2_18[i],jV2_18[i] = nsidc_polar_lonlat(Lon_V2_18[i],Lat_V2_18[i],grid,hemisphere) # V1_17
#
#for i in range(len(Lat_V3_18)):
#    iV3_18[i],jV3_18[i] = nsidc_polar_lonlat(Lon_V3_18[i],Lat_V3_18[i],grid,hemisphere) # V1_17

for i in range(len(Lat_SIPEXII)):
    iSIPEXII[i],jSIPEXII[i] = nsidc_polar_lonlat(Lon_SIPEXII[i],Lat_SIPEXII[i],grid,hemisphere) # V1_17
    
#------------------------------------------------------------------------------
# Save the variables as a .csv

#dfV1_17 = np.column_stack((Date_V1_17,iV1_17,jV1_17,Lon_V1_17,Lat_V1_17))
#dfV1_17 = pd.DataFrame.from_dict(dfV1_17)
#dfV1_17.columns = ['Date_V1_17','I','J','Lon_V1_17','Lat_V1_17']
#dfV1_17.to_csv('/Users/ncp532/Documents/Data/V1_17_M_ij.csv')
#
#dfV2_17 = np.column_stack((Date_V2_17,iV2_17,jV2_17,Lon_V2_17,Lat_V2_17))
#dfV2_17 = pd.DataFrame.from_dict(dfV2_17)
#dfV2_17.columns = ['Date_V2_17','I','J','Lon_V2_17','Lat_V2_17']
#dfV2_17.to_csv('/Users/ncp532/Documents/Data/V2_17_M_ij.csv')
#
#dfV3_17 = np.column_stack((Date_V3_17,iV3_17,jV3_17,Lon_V3_17,Lat_V3_17))
#dfV3_17 = pd.DataFrame.from_dict(dfV3_17)
#dfV3_17.columns = ['Date_V3_17','I','J','Lon_V3_17','Lat_V3_17']
#dfV3_17.to_csv('/Users/ncp532/Documents/Data/V3_17_M_ij.csv')
#
#dfV1_18 = np.column_stack((Date_V1_18,iV1_18,jV1_18,Lon_V1_18,Lat_V1_18))
#dfV1_18 = pd.DataFrame.from_dict(dfV1_18)
#dfV1_18.columns = ['Date_V1_18','I','J','Lon_V1_18','Lat_V1_18']
#dfV1_18.to_csv('/Users/ncp532/Documents/Data/V1_18_M_ij.csv')
#
#dfV2_18 = np.column_stack((Date_V2_18,iV2_18,jV2_18,Lon_V2_18,Lat_V2_18))
#dfV2_18 = pd.DataFrame.from_dict(dfV2_18)
#dfV2_18.columns = ['Date_V2_18','I','J','Lon_V2_18','Lat_V2_18']
#dfV2_18.to_csv('/Users/ncp532/Documents/Data/V2_18_M_ij.csv')
#
#df6 = np.column_stack((Date_V3_18,iV3_18,jV3_18,Lon_V3_18,Lat_V3_18))
#df6 = pd.DataFrame.from_dict(df6)
#df6.columns = ['Date_V3_18','I','J','Lon_V3_18','Lat_V3_18']
#df6.to_csv('/Users/ncp532/Documents/Data/V3_18_M_ij.csv')

df7 = np.column_stack((Date_SIPEXII,iSIPEXII,jSIPEXII,Lon_SIPEXII,Lat_SIPEXII))
df7 = pd.DataFrame.from_dict(df7)
df7.columns = ['Date_SIPEXII','I','J','Lon_SIPEXII','Lat_SIPEXII']
df7.to_csv('/Users/ncp532/Documents/Data/SIPEXII_M_ij.csv')
