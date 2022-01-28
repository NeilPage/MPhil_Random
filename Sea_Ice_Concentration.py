#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 08:32:21 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
import xarray as xr

# DATA HANDLING PACKAGES
import numpy as np
import matplotlib.pyplot as plt

# DRAWING PACKAGES
import cartopy.crs as ccrs
from matplotlib import cm                   # imports the colormap function from matplotlib
import cartopy
#------------------------------------------------------------------------------
# DEFINE THE DATASET

# NETCDF
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20190101_median5day.nc')
#cubes2 = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_icdr_sh_f18_20180101_v01r00.nc')
#cubes2 = MFDataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_icdr_sh_f18_2018*.nc')

# XARRAY
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20190101_median5day.nc',decode_cf=False,engine='netcdf4')
#cubes2 = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_icdr_sh_f18_20180201_v01r00.nc',decode_cf=False,engine='netcdf4')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

# NETCDF
#-----------------------
# Hamburg
#lats = cubes.variables['latitude'][:] # (y,x) (664,632)
#lons = cubes.variables['longitude'][:] # (y,x) (664,632)
#land = cubes.variables['land'][0,:,:] # (time,y,x) (1,664,632)
#seaice_data = cubes.variables['sea_ice_area_fraction'][0,:,:] # Sea Ice concentration (time,y,x)(1,664,632)

# NSIDC
#seaice_data = cubes2.variables['seaice_conc_cdr'][0,:,:]
#lats = cubes2.variables['latitude'][:]
#lons = cubes2.variables['longitude'][:]

# XARRAY
#-----------------------
# Hamburg
#lats = cubes.latitude]
#lons = cubes.longitude
#seaice_data = cubes.sea_ice_area_fraction[0,:,:]

# NSIDC
# seaice_data = cubes2.seaice_conc_cdr[0,:,:]
# lats = cubes2.latitude
# lons = cubes2.longitude

#------------------------------------------------------------------------------
# START PLOTTING THE MAP

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

# We will view the map on a polar stereo coordinate system
projection = ccrs.SouthPolarStereo(central_longitude=0)

plt.close()
#fig = plt.figure(figsize=(6, 3))
fig = plt.figure()

#------------------------------------
# PLOT AX1 (SEA ICE CONCENTRATION)
ax1 = plt.subplot(111, projection=projection)

#ax.set_global()
ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='white')
ax1.coastlines()
gl = ax1.gridlines(color='gray',alpha=0.5,)
ax1.set_extent([-180, 180, -55, -90], crs=ccrs.PlateCarree())

#------------------------------------
# Plot the data 
#cs = ax1.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cm.get_cmap('jet',20)) #, bins=np.arange(0,100, 10))
#cb = fig.colorbar(cs,ticks=[0,10,20,30,40,50,60,70,80,90,100])

# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276

# Plot the station markers
ax1.plot(Davis_lat, Davis_lon, transform=data_crs, color='k', marker='o')
ax1.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='k', marker='o')
ax1.plot(Casey_lat, Casey_lon, transform=data_crs, color='k', marker='o')

# Plot the mareker labels
ax1.text(Davis_lat + 2, Davis_lon - 1, 'Davis', transform=data_crs, horizontalalignment='right')
ax1.text(Mawson_lat + 2, Mawson_lon - 1, 'Mawson', transform=data_crs, horizontalalignment='right')
ax1.text(Casey_lat + 2, Casey_lon - 1, 'Casey', transform=data_crs, horizontalalignment='right')


#------------------------------------
# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("Sea Ice Cover (01/01/2019)")
#cb.set_label('Concentration (%)', rotation=90)
