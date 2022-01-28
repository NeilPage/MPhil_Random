#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:22:07 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
#from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
from netCDF4 import Dataset,MFDataset
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
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20180112_median5day.nc')

# NSIDC
# 2017
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20180103_v03r01.nc')
# 2018
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/seaice_conc_daily_icdr_sh_f18_20181230_v01r00.nc')

# XARRAY
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20180215_median5day.nc',decode_cf=False,engine='netcdf4')
# 2017-18
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20180201_v03r01.nc',decode_cf=False,engine='netcdf4')
# 2018-19
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_icdr_sh_f18_20180201_v01r00.nc',decode_cf=False,engine='netcdf4')
# SIPEXII
cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2012/seaice_conc_daily_sh_f17_20120912_v03r01.nc',decode_cf=False,engine='netcdf4')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

# NETCDF
#lats = cubes.variables['latitude'][:] # (y,x) (664,632)
#lons = cubes.variables['longitude'][:] # (y,x) (664,632)
# NSIDC
#seaice_data = cubes.variables['seaice_conc_cdr'][0,:,:]*100 # Sea Ice concentration (time,y,x)(1,664,632)
# HAMBURG
#seaice_data = cubes.variables['sea_ice_area_fraction'][0,:,:] # Sea Ice concentration (time,y,x)(1,664,632)

# XARRAY
lats = cubes.latitude
lons = cubes.longitude
#land = cubes.land
#seaice_data = cubes.sea_ice_area_fraction[0,:,:] # Hamburg
seaice_data = cubes.seaice_conc_cdr[0,:,:] # Hamburg

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

#------------------------------------------------------------------------------
# START PLOTTING THE MAP
fig = plt.figure(figsize=(6, 3))
plt.subplots_adjust()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([50, 70, -60, -70])#, crs=ccrs.PlateCarree())
#ax2 = ax.twinx()
#ax2.set_extent([40, 120, -50, -75])#, crs=ccrs.PlateCarree())

ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='white')
ax.coastlines()

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
gl.xlocator = mticker.FixedLocator([40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,110,115,120,125,130,135,140,145,150,155,160])
gl.ylocator = mticker.FixedLocator([-55,-57.5,-60,-62.5,-65,-67.5,-70,-72.5,-75,-77.5,-80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#------------------------------------
# PLOT THE DATA (SEA ICE CONCENTRATION) 

cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cm.get_cmap('jet',20)) #, bins=np.arange(0,100, 10))
cb = fig.colorbar(cs,ticks=[0,10,20,30,40,50,60,70,80,90,100],shrink=.80)
#cb = fig.colorbar(cs,ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],shrink=.50)
#------------------------------------
# PLOT THE AAD STATIONS

# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276
#Event1_lon, Event1_lat = -66.2699, 110.5394
#Event2_lon, Event2_lat = -68.5741, 77.9304
Event3_lon, Event3_lat = -67.3525, 62.8075
#Event4_lon, Event4_lat = -68.5650, 77.8442
#Event5_lon, Event5_lat = -63.7395, 81.6975
#Event6_lon, Event6_lat = -67.5688, 62.7511
#Event7_lon, Event7_lat = -64.3930, 120.0869

# Plot the station markers
ax.plot(Davis_lat, Davis_lon, color='k', marker='o')
ax.plot(Mawson_lat, Mawson_lon, color='k', marker='o')
ax.plot(Casey_lat, Casey_lon, color='k', marker='o')
#ax.plot(Event1_lat, Event1_lon, color='k', marker='o')
#ax.plot(Event2_lat, Event2_lon, color='k', marker='o')
ax.plot(Event3_lat, Event3_lon, color='k', marker='o')
#ax.plot(Event4_lat, Event4_lon, color='k', marker='o')
#ax.plot(Event5_lat, Event5_lon, color='k', marker='o')
#ax.plot(Event6_lat, Event6_lon, color='k', marker='o')
#ax.plot(Event7_lat, Event7_lon, color='k', marker='o')

# Plot the mareker labels
ax.text(Davis_lat + 1.5, Davis_lon - 0.5, 'Davis',horizontalalignment='right', fontsize=15)
ax.text(Mawson_lat + 1.5, Mawson_lon - 0.5, 'Mawson',horizontalalignment='right', fontsize=15)
ax.text(Casey_lat + 1.5, Casey_lon - 0.5, 'Casey',horizontalalignment='right', fontsize=15)
#ax.text(Event1_lat + 1.5, Event1_lon + 0.5, 'V2 (2017-18)',horizontalalignment='right', fontsize=15)
#ax.text(Event2_lat + 1.5, Event2_lon + 0.5, 'V3 (2017-18)',horizontalalignment='right', fontsize=15)
ax.text(Event3_lat + 1.5, Event3_lon + 0.5, 'V3 (2017-18)',horizontalalignment='right', fontsize=15)
#ax.text(Event4_lat + 1.5, Event4_lon + 0.5, 'V1 (2018-19)',horizontalalignment='right', fontsize=15)
#ax.text(Event5_lat + 1.5, Event5_lon - 0.5, 'V3 (2018-19)',horizontalalignment='right', fontsize=15)
#ax.text(Event6_lat + 1.5, Event6_lon + 0.5, 'V3 (2018-19)',horizontalalignment='right', fontsize=15)
#ax.text(Event7_lat + 1.5, Event7_lon - 0.5, 'SIPEXII',horizontalalignment='right', fontsize=15)

#------------------------------------
# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("Sea Ice Cover (1 Febuary 2018)", y=1.1, fontsize=20)
cb.set_label('Concentration (%)', rotation=90, fontsize=15)

ax.text(-0.04, 0.55, 'latitude [$^\circ$]', fontsize=15, va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)

ax.text(0.5, -0.1, 'longitude [$^\circ$]', fontsize=15, va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
plt.show()