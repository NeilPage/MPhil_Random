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

# NETCDF
cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20171231_median5day.nc')
#seaice = cubes.variables['sea_ice_area_fraction'][0,:,:]

# XARRAY
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20190101_median5day.nc',decode_cf=False,engine='netcdf4')
#seaice = cubes.sea_ice_area_fraction[0,:,:]

# NSIDC
#cubes2 = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_icdr_sh_f18_20180201_v01r00.nc',decode_cf=False,engine='netcdf4')
#cubes2 = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/seaice_conc_daily_icdr_sh_f18_20181120_v01r00.nc')
#cubes2 = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20171231_v03r01.nc')

V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V1_17_Hourly_LatLong.csv')
V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V2_17_Hourly_LatLong.csv')
V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/V3_17_Hourly_LatLong.csv')
V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1_18_Hourly_LatLong.csv')
V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V2_18_Hourly_LatLong.csv')
V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V3_18_Hourly_LatLong.csv')

Date_V1_17 = np.array(V1_17['Date_V1_17'])
Lat_V1_17 = np.array(V1_17['Lat_V1_17'])
Lon_V1_17 = np.array(V1_17['Lon_V1_17'])
Date_V2_17 = np.array(V2_17['Date_V2_17'])
Lat_V2_17 = np.array(V2_17['Lat_V2_17'])
Lon_V2_17 = np.array(V2_17['Lon_V2_17'])
Date_V3_17 = np.array(V3_17['Date_V3_17'])
Lat_V3_17 = np.array(V3_17['Lat_V3_17'])
Lon_V3_17 = np.array(V3_17['Lon_V3_17'])
Date_V1_18 = np.array(V1_18['Date_V1_18'])
Lat_V1_18 = np.array(V1_18['Lat_V1_18'])
Lon_V1_18 = np.array(V1_18['Lon_V1_18'])
Date_V2_18 = np.array(V2_18['Date_V2_18'])
Lat_V2_18 = np.array(V2_18['Lat_V2_18'])
Lon_V2_18 = np.array(V2_18['Lon_V2_18'])
Date_V3_18 = np.array(V3_18['Date_V3_18'])
Lat_V3_18 = np.array(V3_18['Lat_V3_18'])
Lon_V3_18 = np.array(V3_18['Lon_V3_18'])

#------------------------------------------------------------------------------
# DEFINE THE NSIDC VARIABLES

# NETCDF
#seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2018-19
#seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:] # NSIDC 2017
seaice_data = cubes.variables['sea_ice_area_fraction'][0,:,:] # Hamburg
y = cubes.variables['latitude'][:] # (y,x) (664,632)
x = cubes.variables['longitude'][:] # (y,x) (664,632)

# look up nearest grid square
#davis = cubes2.sel(0,x=-68.5766,y=77.9674,method='nearest')
a = seaice_data[512][422] 

# XARRAY
#seaice_data = cubes2.seaice_conc_cdr[0,:,:] # y = 332, x = 316
#y = cubes2.latitude
#x = cubes2.longitude

true_scale_lat = 70 # true scale latitude in degrees (NSIDC made this polar stereo projection true at 70 degrees rather than at the poles)
re = 6378.273 # Earth radius in km
e =  0.081816153 # Earth eccentricity
hemisphere = -1 # 1 for Northern hemisphere, -1 for Southern
grid = 25

#------------------------------------------------------------------------------
# FUNCTION TO CONVERT FROM POLAR STEREOGRAPHIC (x,y) TO LAT/LONG COORDINATES

def polar_xy_to_lonlat(x, y, true_scale_lat, re, e, hemisphere):
    """Convert from Polar Stereographic (x, y) coordinates to
    geodetic longitude and latitude.
    Args:
        x (float): X coordinate(s) in km
        y (float): Y coordinate(s) in km
        true_scale_lat (float): true-scale latitude in degrees
        hemisphere (1 or -1): 1 for Northern hemisphere, -1 for Southern
        re (float): Earth radius in km
        e (float): Earth eccentricity
    Returns:
        If x and y are scalars then the result is a
        two-element list containing [longitude, latitude].
        If x and y are numpy arrays then the result will be a two-element
        list where the first element is a numpy array containing
        the longitudes and the second element is a numpy array containing
        the latitudes.
    """

    e2 = e * e
    slat = true_scale_lat * np.pi / 180
    rho = np.sqrt(x ** 2 + y ** 2)

    if abs(true_scale_lat - 90.) < 1e-5:
        t = rho * np.sqrt((1 + e) ** (1 + e) * (1 - e) ** (1 - e)) / (2 * re)
    else:
        cm = np.cos(slat) / np.sqrt(1 - e2 * (np.sin(slat) ** 2))
        t = np.tan((np.pi / 4) - (slat / 2)) / \
            ((1 - e * np.sin(slat)) / (1 + e * np.sin(slat))) ** (e / 2)
        t = rho * t / (re * cm)

    chi = (np.pi / 2) - 2 * np.arctan(t)
    lat = chi + \
        ((e2 / 2) + (5 * e2 ** 2 / 24) + (e2 ** 3 / 12)) * np.sin(2 * chi) + \
        ((7 * e2 ** 2 / 48) + (29 * e2 ** 3 / 240)) * np.sin(4 * chi) + \
        (7 * e2 ** 3 / 120) * np.sin(6 * chi)
    lat = hemisphere * lat * 180 / np.pi
    lon = np.arctan2(hemisphere * x, -hemisphere * y)
    lon = hemisphere * lon * 180 / np.pi
    lon = lon + np.less(lon, 0) * 360
    return [lon, lat]

#------------------------------------------------------------------------------
# FUNCTION TO CONVERT FROM LAT/LONG COORDINATES TO POLAR STEREOGRAPHIC (x,y)
    
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

#------------------------------------------------------------------------------
# FUNCTION TO TRANSFORM FROM NSIDC POLAR STEREOGRAPHIC (I,J) TO LAT/LONG COORDINATES
    
def nsidc_polar_ij(i, j, grid, hemisphere):
    """Transform from NSIDC Polar Stereographic I, J coordinates
    to longitude and latitude coordinates
    
    Args:
        i (int): an integer or integer array giving the x grid coordinate(s)
        j (int): an integer or integer array giving the y grid coordinate(s)
        grid (float): 6.25, 12.5 or 25; the grid cell dimensions in km
        hemisphere (1 or -1): Northern or Southern hemisphere
    
    Returns:
        If i and j are scalars then the result is a
        two-element list containing [longitude, latitude].
        If i and j are numpy arrays then the result will be a two-element
        list where the first element is a numpy array containing
        the longitudes and the second element is a numpy array containing
        the latitudes.
    Examples:
        print(nsidc_polar_ij(608, 896, 12.5, 1))
            [350.01450147320855, 34.40871032516291]
    """

    true_scale_lat = 70
    re = 6378.273
    e = 0.081816153

    if grid != 6.25 and grid != 12.5 and grid != 25:
        raise ValueError("Legal grid values are 6.25, 12.5, or 25")
    
    if hemisphere != 1 and hemisphere != -1:
        raise ValueError("Legal hemisphere values are 1 or -1")

    if hemisphere == 1:
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

    if np.any(np.less(i, 1)) or np.any(np.greater(i, imax)):
        raise ValueError("'i' value is out of range: [1, " + str(imax) + "]")
    if np.any(np.less(j, 1)) or np.any(np.greater(j, jmax)):
        raise ValueError("'j' value is out of range: [1, " + str(jmax) + "]")

    # Convert I, J pairs to x and y distances from origin.
    x = ((i - 1)*grid) + xmin
    y = ((jmax - j)*grid) + ymin
    lonlat = polar_xy_to_lonlat(x, y, true_scale_lat, re, e, hemisphere)
    lon = lonlat[0] - delta
    lon = lon + np.less(lon, 0)*360
    return [lon, lonlat[1]]

#------------------------------------------------------------------------------
# FUNCTION TO TRANSFORM FROM LAT/LONG COORDINATES TO NSIDC POLAR STEREOGRAPHIC (I,J)
  
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

# polar_xy_to_lonlat
lon,lat = polar_xy_to_lonlat(x, y, true_scale_lat, re, e, hemisphere) 

i,j = nsidc_polar_lonlat(77.9674, -68.5766, grid, hemisphere) # Davis
i,j = nsidc_polar_lonlat(110.5276, -66.2818, grid, hemisphere) # Casey
i,j = nsidc_polar_lonlat(62.8738, -67.6027, grid, hemisphere) # Mawson

 # Set up reservoir arrays
iV1_17 = np.zeros(len(Lat_V1_17))
jV1_17 = np.zeros(len(Lat_V1_17))
iV2_17 = np.zeros(len(Lat_V2_17))
jV2_17 = np.zeros(len(Lat_V2_17))
iV3_17 = np.zeros(len(Lat_V3_17))
jV3_17 = np.zeros(len(Lat_V3_17))
iV1_18 = np.zeros(len(Lat_V1_18))
jV1_18 = np.zeros(len(Lat_V1_18))
iV2_18 = np.zeros(len(Lat_V2_18))
jV2_18 = np.zeros(len(Lat_V2_18))
iV3_18 = np.zeros(len(Lat_V3_18))
jV3_18 = np.zeros(len(Lat_V3_18))

for i in range(len(Lat_V1_17)):
    iV1_17[i],jV1_17[i] = nsidc_polar_lonlat(Lon_V1_17[i],Lat_V1_17[i],grid,hemisphere) # V1_17

for i in range(len(Lat_V2_17)):
    iV2_17[i],jV2_17[i] = nsidc_polar_lonlat(Lon_V2_17[i],Lat_V2_17[i],grid,hemisphere) # V2_17

for i in range(len(Lat_V3_17)):
    iV3_17[i],jV3_17[i] = nsidc_polar_lonlat(Lon_V3_17[i],Lat_V3_17[i],grid,hemisphere) # V1_17

for i in range(len(Lat_V1_18)):
    iV1_18[i],jV1_18[i] = nsidc_polar_lonlat(Lon_V1_18[i],Lat_V1_18[i],grid,hemisphere) # V1_17

for i in range(len(Lat_V2_18)):
    iV2_18[i],jV2_18[i] = nsidc_polar_lonlat(Lon_V2_18[i],Lat_V2_18[i],grid,hemisphere) # V1_17

for i in range(len(Lat_V3_18)):
    iV3_18[i],jV3_18[i] = nsidc_polar_lonlat(Lon_V3_18[i],Lat_V3_18[i],grid,hemisphere) # V1_17

dfV1_17 = np.column_stack((Date_V1_17,iV1_17,jV1_17,Lon_V1_17,Lat_V1_17))
dfV1_17 = pd.DataFrame.from_dict(dfV1_17)
dfV1_17.columns = ['Date_V1_17','I','J','Lon_V1_17','Lat_V1_17']
dfV1_17.to_csv('/Users/ncp532/Documents/Data/V1_17_H_ij.csv')

dfV2_17 = np.column_stack((Date_V2_17,iV2_17,jV2_17,Lon_V2_17,Lat_V2_17))
dfV2_17 = pd.DataFrame.from_dict(dfV2_17)
dfV2_17.columns = ['Date_V2_17','I','J','Lon_V2_17','Lat_V2_17']
dfV2_17.to_csv('/Users/ncp532/Documents/Data/V2_17_H_ij.csv')

dfV3_17 = np.column_stack((Date_V3_17,iV3_17,jV3_17,Lon_V3_17,Lat_V3_17))
dfV3_17 = pd.DataFrame.from_dict(dfV3_17)
dfV3_17.columns = ['Date_V3_17','I','J','Lon_V3_17','Lat_V3_17']
dfV3_17.to_csv('/Users/ncp532/Documents/Data/V3_17_H_ij.csv')

dfV1_18 = np.column_stack((Date_V1_18,iV1_18,jV1_18,Lon_V1_18,Lat_V1_18))
dfV1_18 = pd.DataFrame.from_dict(dfV1_18)
dfV1_18.columns = ['Date_V1_18','I','J','Lon_V1_18','Lat_V1_18']
dfV1_18.to_csv('/Users/ncp532/Documents/Data/V1_18_H_ij.csv')

dfV2_18 = np.column_stack((Date_V2_18,iV2_18,jV2_18,Lon_V2_18,Lat_V2_18))
dfV2_18 = pd.DataFrame.from_dict(dfV2_18)
dfV2_18.columns = ['Date_V2_18','I','J','Lon_V2_18','Lat_V2_18']
dfV2_18.to_csv('/Users/ncp532/Documents/Data/V2_18_H_ij.csv')

df6 = np.column_stack((Date_V3_18,iV3_18,jV3_18,Lon_V3_18,Lat_V3_18))
df6 = pd.DataFrame.from_dict(df6)
df6.columns = ['Date_V3_18','I','J','Lon_V3_18','Lat_V3_18']
df6.to_csv('/Users/ncp532/Documents/Data/V3_18_H_ij.csv')

lon,lat = nsidc_polar_ij(i, j, grid, hemisphere)

# nsidc_polar_ij
#lon,lonlat = nsidc_polar_ij(x, y, grid, hemisphere)

#lon = lon[:,0] # longitude data
#lat = lat[:,0] # lattude data

lonlat = np.column_stack((lon,lat))
lonlat = pd.DataFrame.from_dict(lonlat)

#------------------------------------------------------------------------------
# GET THE DATA FROM THE GRID SQUARE NEAREST TO EACH AAD STATION

# Coordinates of AAD stations (lat/lon) and (WGS84 Antarctic Polar Stereographic (EPSG:3031))
davis_lat, davis_lon = -68.5766, 77.9674 # Y (490757.62762) X (2302390.776809)
casey_lon, casey_lat = -66.2818, 110.5276 # Y (-916269.355788) X (2447079.532811)
mawson_lon, mawson_lat = -67.6027, 62.8738 # Y (1123346.333303) X (2192738.114197)

# look up nearest grid square
#davis = cubes2.sel(0,x=-68.5766,y=77.9674,method='nearest')
#a = seaice_data[250][155] # Davis
#b = seaice_data[256][211] # Casey
#c = seaice_data[246][130] # Mawson