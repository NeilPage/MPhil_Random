#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:22:07 2019

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
import matplotlib.ticker as mticker
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# HAMBURG
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20181215_median5day.nc')

# NSIDC
# 2017
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20180105_v03r01.nc')
# 2018
cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/seaice_conc_daily_icdr_sh_f18_20181230_v01r00.nc')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

lats = cubes.variables['latitude'][:] # (y,x) (664,632)
lons = cubes.variables['longitude'][:] # (y,x) (664,632)
# NSIDC
seaice_data = cubes.variables['seaice_conc_cdr'][0,:,:] # Sea Ice concentration (time,y,x)(1,664,632)
# HAMBURG
#seaice_data = cubes.variables['sea_ice_area_fraction'][0,:,:] # Sea Ice concentration (time,y,x)(1,664,632)

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

#------------------------------------------------------------------------------
# START PLOTTING THE MAP
fig = plt.figure(figsize=(6, 3))
plt.subplots_adjust()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([100, 120, -57.5, -72.5])#, crs=ccrs.PlateCarree())
#ax2 = ax.twinx()
#ax2.set_extent([40, 120, -50, -75])#, crs=ccrs.PlateCarree())

ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='white')
ax.coastlines()

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
gl.xlocator = mticker.FixedLocator([100,104,108,112,116,120])
gl.ylocator = mticker.FixedLocator([-57.5,-60,-62.5,-65,-67.5,-70,-72.5,-75])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#------------------------------------
# PLOT THE DATA (SEA ICE CONCENTRATION) 

cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cm.get_cmap('jet',20)) #, bins=np.arange(0,100, 10))
cb = fig.colorbar(cs,ticks=[0,10,20,30,40,50,60,70,80,90,100])

#------------------------------------
# PLOT THE AAD STATIONS

# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276 # (-66.2818,110.5276)

# Plot the station markers
ax.plot(Davis_lat, Davis_lon, color='k', marker='o')
ax.plot(Mawson_lat, Mawson_lon, color='k', marker='o')
ax.plot(Casey_lat, Casey_lon, color='k', marker='o')

# Plot the mareker labels
ax.text(Davis_lat + 3, Davis_lon - 2, 'Davis',horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson',horizontalalignment='right')
ax.text(Casey_lat + 1, Casey_lon - 1, 'Casey',horizontalalignment='right')

#------------------------------------
# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("Sea Ice Cover (30/12/2018)", y=1.05, fontsize=20)
cb.set_label('Concentration (%)', rotation=90)

ax.text(-0.065, 0.55, 'latitude [$^\circ$]', fontsize=15, va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)

ax.text(0.5, -0.1, 'longitude [$^\circ$]', fontsize=15, va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
plt.show()