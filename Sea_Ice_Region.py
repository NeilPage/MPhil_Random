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
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt

# DRAWING PACKAGES
import cartopy.crs as ccrs
from matplotlib import cm, transforms                  # imports the colormap function from matplotlib
import matplotlib.ticker as mticker
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker

#------------------------------------------------------------------------------
# DEFINE THE DATASET

#-------------
# Sea Ice
#-------------

# XARRAY
cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20181107_median5day.nc',decode_cf=False,engine='netcdf4')
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20181215_median5day.nc',decode_cf=False,engine='netcdf4')
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2012/seaice_conc_daily_sh_f17_20121011_v03r01.nc',decode_cf=False,engine='netcdf4')

#-------------
# Back Trajectories
#-------------
#Traj = pd.read_csv('/Users/ncp532/Documents/data/SeaIce_Trajectories/gdas1nov2000spring2018110701.csv',index_col=0)

# Set the location for the working directory
os.chdir("/Users/ncp532/Documents/Data/SeaIce_Trajectories/100m/")

# Set the file type to .csv
extension = 'csv'
# Save a list of all the file names to the variable all_filenames 
all_filenames = [i for i in glob.glob('gdas1nov0100spring20181107*.{}'.format(extension))]
#all_filenames = [i for i in glob.glob('gdas1dec0100summer20181215*.{}'.format(extension))]

# Combine all files in the list
Traj = pd.concat([pd.read_csv(f) for f in all_filenames ])


# # Seperate Dataframes
# d = {}  # dictionary that will hold them 
# for f in all_filenames:  # loop over files
#    # read csv into a dataframe and add it to dict with file_name as it key
#    d[f] = pd.read_csv(f)


#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

# NETCDF
#lats = cubes.variables['latitude'][:] # (y,x) (664,632)
#lons = cubes.variables['longitude'][:] # (y,x) (664,632)
# NSIDC
#seaice_data = cubes.variables['seaice_conc_cdr'][0,:,:]*100 # Sea Ice concentration (time,y,x)(1,664,632)
# HAMBURG
seaice_data = cubes.variables['sea_ice_area_fraction'][0,:,:] # Sea Ice concentration (time,y,x)(1,664,632)
mask        = cubes.variables['land'][0,:,:]                  # land mask
#mask2 = np.ma.masked_less(mask,1)
mask2 = np.ma.masked_where(mask == 1, mask)

# XARRAY
lats = cubes.latitude
lons = cubes.longitude
#land = cubes.land
#seaice_data = cubes.sea_ice_area_fraction[0,:,:] # Hamburg
#seaice_data = cubes.seaice_conc_cdr[0,:,:] # Hamburg

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

#------------------------------------------------------------------------------
# PLOT THE FIGURE

fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
# gs = gridspec.GridSpec(nrows=3,
#                        ncols=4, 
#                        figure=fig, 
#                        width_ratios= [0.25,0.25,0.25,0.25],
#                        height_ratios=[0.25, 0.25, 0.25],
#                        hspace=0.3, wspace=0.35)

#-------------------------------------
# SUBPLOT 1 (MAP)
#ax = plt.subplot(gs[0,0:2])
ax = plt.subplot(211,projection=ccrs.PlateCarree()) # options graph 1 (vertical no, horizontal no, graph no)
#ax = plt.subplot(211,projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no, graph no)

#ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([40, 140, -75, -55])#, crs=ccrs.PlateCarree())
#ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())

ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

#------------------------------------
# PLOT THE DATA (SEA ICE CONCENTRATION) 

cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cm.get_cmap('viridis',20)) #, bins=np.arange(0,100, 10))

cmap1 = cm.get_cmap("Greys",lut=10)
cmap1.set_bad("grey")

#ax.pcolormesh(lons, lats, mask2, transform=data_crs, cmap=cmap1, alpha=0)
cb = fig.colorbar(cs,ticks=[0,10,20,30,40,50,60,70,80,90,100])#,shrink=.50)
#cb = plt.colorbar(cs,ticks=[0,10,20,30,40,50,60,70,80,90,100])

#------------------------------------
# PLOT THE BACK TRAJECTORIES

# Traj 1
ax.plot(Traj['Traj Lon'],Traj['Traj Lat'], c='magenta', transform=data_crs, linewidth=0.5)#, s=10, label="SIPEXII") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  

#------------------------------------
# PLOT THE AAD STATIONS

# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat, Davis_lon, transform=data_crs, color='k', marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='k', marker='*')
ax.plot(Casey_lat, Casey_lon, transform=data_crs, color='k', marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3, Davis_lon - 2, 'Davis', transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3, Casey_lon - 2, 'Casey',transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')


gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
gl.xlocator = mticker.FixedLocator([-180,-170,-160,-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
gl.ylocator = mticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#------------------------------------
# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (07/11/2018)", y=1.1, fontsize=20)
plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)
cb.set_label('Concentration (%)', rotation=90)

ax.text(-0.07, 0.55, 'latitude [$^\circ$]', fontsize=15, va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)

ax.text(0.5, -0.2, 'longitude [$^\circ$]', fontsize=15, va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)

#-------------------------------------
# SUBPLOT 1 (MAP)
#ax = plt.subplot(gs[1,:2])
ax = plt.subplot(212) # options graph 1 (vertical no, horizontal no, graph no)

# Plot axis lines for height
plt.axhline(500,  linewidth=0.5, color='k')
plt.axhline(1000, linewidth=0.5, color='k')
plt.axhline(1500, linewidth=0.5, color='k')
plt.axhline(2000, linewidth=0.5, color='k')
plt.axhline(2500, linewidth=0.5, color='k')
plt.axhline(3000, linewidth=0.5, color='k')

# Back trajectory altitude
ax.plot(Traj['Traj Age'], Traj['Traj Height (m)'], marker='o', c='r', markersize = 2.0, linestyle='none', label='Traj Height (m)')
ax.plot(Traj['Traj Age'][0], Traj['Traj Height (m)'][0], color='k', marker='*')

# Format x-axis
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
#ax.set_xlim(,0)

# Format y-axis
ax.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))
ax.set_ylim(0,)

# Plot axis labels & title
plt.title("Back trajectory height", y=1.1, fontsize=20)
ax.set_xlabel('Age (hours)', fontsize=15)
ax.set_ylabel('Height (MSL)', fontsize=15)
#cb = fig.colorbar(cs, ax=ax, extend='max', pad = 0.16).ax.set_visible(False)

