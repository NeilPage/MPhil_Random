#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 11:56:27 2020

@author: ncp532
"""

# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files
import xarray as xr

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats, interpolate

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# DEFINE THE DATASETS

#--------------
# MERRA-2
#--------------

# A1
A1_V1_17   = xr.open_mfdataset('/Users/ncp532/Documents/Data/MERRA2/V3_18/MERRA2.201*.A1.2x25.nc4', combine='by_coords') # A1 Fields

# A3
A3_V1_17   = xr.open_mfdataset('/Users/ncp532/Documents/Data/MERRA2/V3_18/MERRA2.201*.A3dyn.2x25.nc4', combine='by_coords') # A3dyn Fields

# I3
I3_V1_17   = xr.open_mfdataset('/Users/ncp532/Documents/Data/MERRA2/V3_18/MERRA2.201*.I3.2x25.nc4', combine='by_coords') # I3 Fields

#------------------------------------------------------------------------------
# RESAMPLE THE I3 & A3 DATASETS FROM 3 HOURS TO 1 HOUR, A1 FROM MIDPOINT TO START

A1_V1_17 = A1_V1_17.resample(time="60T").mean()

# A3
A3_V1_17 = A3_V1_17.resample(time="60T").mean()

# I3
I3_V1_17 = I3_V1_17.resample(time="60T").mean()

#------------------------------------------------------------------------------
# INTERPOLATE THE VALUES WITHIN THE DATASETS
#----------------------------------
# V1_17 Davis (Lat = 11, Lon = 103)
#----------------------------------

#I3
Temp_V1_17  = pd.DataFrame(I3_V1_17.T[:,:,12,103].values, index=I3_V1_17.time.values).interpolate(method='time')
PS_V1_17    = pd.DataFrame(I3_V1_17.PS[:,12,103].values,  index=I3_V1_17.time.values).interpolate(method='time')

# A3dyn
RH_V1_17    = pd.DataFrame(A3_V1_17.RH[:,:,12,103].values, index=A3_V1_17.time.values).interpolate(method='time')
EWind_V1_17 = pd.DataFrame(A3_V1_17.U[:,:,12,103].values,  index=A3_V1_17.time.values).interpolate(method='time')
NWind_V1_17 = pd.DataFrame(A3_V1_17.V[:,:,12,103].values,  index=A3_V1_17.time.values).interpolate(method='time')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES
#----------------------------------
# V1_17 Davis (Lat = 11, Lon = 103)
#----------------------------------
# A1
Time_V1_17       = A1_V1_17.time[1:1245].values           # DateTime (YYYY-MM_DD HH:MM:SS)
PrecTot_V1_17    = A1_V1_17.PRECTOT[1:1245,12,103].values # Total precipitation (kg/m2/s)
PrecSno_V1_17    = A1_V1_17.PRECSNO[1:1245,12,103].values # Snowfall (kg/m2/s)
Rainfall_V1_17   = PrecTot_V1_17 - PrecSno_V1_17          # Rainfall (kg/m2/s)
SLP_V1_17        = A1_V1_17.SLP[1:1245,12,103].values     # Sea level pressure (Pa)
Temp2m_V1_17     = A1_V1_17.T2M[1:1245,12,103].values     # 2m air temperature (K)
Temp10m_V1_17    = A1_V1_17.T10M[1:1245,12,103].values    # 10m air temperature (K)
QV2m_V1_17       = A1_V1_17.QV2M[1:1245,12,103].values    # 2m specific humidity (kg/kg)
EWind10m_V1_17   = A1_V1_17.U10M[1:1245,12,103].values    # 10m eastward wind (m/s)
NWind10m_V1_17   = A1_V1_17.V10M[1:1245,12,103].values    # 10m northward wind (m/s)
WS10m_V1_17      = np.sqrt(EWind10m_V1_17**2
                           + NWind10m_V1_17**2)           # 10m wind speed (m/s)
# I3
Time2_V1_17      = I3_V1_17.time[1:1245].values           # DateTime (YYYY-MM_DD HH:MM:SS)
Temp_V1_17       = np.array(Temp_V1_17)[1:1245,:]         # Air temperature (K)
PS_V1_17         = np.array(PS_V1_17)[1:1245]             # Surface pressure (Pa)

# A3dyn
Time3_V1_17      = A3_V1_17.time[0:1244].values           # DateTime (YYYY-MM_DD HH:MM:SS)
RH_V1_17         = np.array(RH_V1_17)[0:1244,:]           # Relative humidity (1)
EWind_V1_17      = np.array(EWind_V1_17)[0:1244,:]        # Eastward wind (m/s)
NWind_V1_17      = np.array(NWind_V1_17)[0:1244,:]        # Northward wind (m/s)

##------------------------------------------------------------------------------
## Sanity Plot 1: Time Interpolation
#
#I3Test   = xr.open_mfdataset('/Users/ncp532/Documents/Data/MERRA2/V1_17/MERRA2.201*.I3.2x25.nc4', combine='by_coords') # I3 Fields
#TimeTest = I3Test.time[2:1005].values # DateTime (YYYY-MM_DD HH:MM:SS)
#TempTest = I3Test.T[2:1005,:,79,103].values # DateTime (YYYY-MM_DD HH:MM:SS)
#
#fig = plt.figure()
#ax=plt.subplot(111)
#
## Plot Variables
#ax.plot(TimeTest,    TempTest[:,0], marker='x', c='blue', markersize = 3.0, linestyle='-', label='Temp Original')
#ax.plot(Time2_V1_17, Temp_V1_17[:,0]-5, marker='x', c='red', markersize = 3.0, linestyle='-', label='Temp Interpolated\n(= Original - 5$^\circ$K)')
#
## Format x-axis
#xmajor_formatter = mdates.DateFormatter('%b %Y') # format how the date is displayed
#ax.xaxis.set_major_formatter(xmajor_formatter)
#xminor_formatter = mdates.DateFormatter('%d') # format how the date is displayed
#ax.xaxis.set_minor_formatter(xminor_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#ax.tick_params(axis='x',pad=15)
#
## Format y-axis 1
#ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
#
## Plot the axis labels, legend and title
#plt.title('Sanity Test 1: Time Interpolation', fontsize=25, y=1.02)
#ax.set_ylabel('Temperature (K)', fontsize=15)
#ax.set_xlabel('Date', fontsize=15)
#plt.legend(bbox_to_anchor=(0.85, 0.95), loc=2, borderaxespad=0.)

#------------------------------------------------------------------------------
# CALCULATE THE PRESSURE (GRID CELL CENTRE)

# define Ap (hPa)
Ap = np.array([ 0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
                1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
                4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
                7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
                1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
                1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
                2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
                2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
                1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
                7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01,
                4.017541e+01, 3.381001e+01, 2.836781e+01, 2.373041e+01,
                1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
                9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00,
                4.076571e+00, 3.276431e+00, 2.620211e+00, 2.084970e+00,
                1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
                6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01,
                2.113490e-01, 1.594950e-01, 1.197030e-01, 8.934502e-02,
                6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
                1.000000e-02,]) 

# define Bp (unitless)
Bp = np.array([ 1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
                9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
                8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
                7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
                6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
                4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
                2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
                6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                0.000000e+00,])
        
# variables
#Psurface = PS_V1_17/100 # the "true" surface pressure (Pa/100=hPa) [Time=1003]
Psurface = SLP_V1_17/100 # the "true" surface pressure (Pa/100=hPa) [Time=1003]
Ap       = Ap    # a constant given at level edges (hPa) [Levels=73]
Bp       = Bp    # a constant given at level edges (unitless) [Levels=73]
Ptop     = 0.010 # The pressure at the model top (hPa)

# set an empty array for PedgeL, ETAedge and ETAcentre
PedgeL    = np.zeros([len(Psurface),len(Ap)])   # hPa [Time=1003,Levels=73]
ETAedge   = np.zeros([len(Psurface),len(Ap)])   # [Time=1003,Levels=73]
ETAcentre = np.zeros([len(Psurface),len(Ap)-1]) # [Time=1003,Levels=72]

#----------
# equations
#----------

# calculate PedgeL (Pedge=Ap+[Bp*Psurface])
for i in range(0,len(Psurface),1):
    PedgeL[i]  = Ap + (Bp * Psurface[i]) # pressure at bottom edge of each level (hPa)

# calcuate PedgeL2     
PedgeL2  = PedgeL[:,1:73] # pressure at bottom edge of level above (hPa)
    
# calculate PcentreL (Pcentre=[Pedge+Pedge+1]/2)
PcentreL  = (PedgeL[:,0:72] + PedgeL2)/2 # pressure at the centre of the level (hPa)

# calculate ETAedge (ETAedge=[Pedge-Ptop]/[Psurface-Ptop])
for i in range(0,len(Psurface),1):
    ETAedge[i] = (PedgeL[i] - Ptop) / (Psurface[i] - Ptop) # ETA coordinate edges

# calculate ETAcentre (ETAcentre=[Pcentre-Ptop]/[Psurface-Ptop])
for i in range(0,len(Psurface),1):
    ETAcentre[i] = (PcentreL[i] - Ptop) / (Psurface[i] - Ptop) # ETA coordinate centres

#------------------------------------------------------------------------------
# CALCULATE THE ALTITUDE (GRID CELL CENTRE)

# variables
PcentreL   = PcentreL   # pressure at the centre of the level (hPa)
Psurface   = Psurface   # the "true" surface pressure (hPa)
g          = 9.80665    # acceleration due to gravity (m/s2)
Mda        = 0.02897    # molar mass dry air (kg/mol)
R          = 8.3143     # specific gas constant (J/K/mol)
Temp_V1_17 = Temp_V1_17 # temperature at the centre of the level (K)

# set an empty array for h & h_RH
h    = np.zeros([len(Psurface),len(Ap)-1]) # [Time=1003,Levels=72]
h_RH = np.zeros([len(Psurface),len(Ap)-1]) # [Time=1003,Levels=72] 

#----------
# equations
#----------

# barometric equation (calculate pressure at a given altitude)
# p = p0*exp(-1*((Mdryair*g)/(R*T)*h))

# rearrange the barometric equation to solve for h (altitude at a given pressure)
# h = -1*((ln(p/p0)*R*T)/(Mdryair*g))

# equation (np.log is the same as ln)
for i in range(0,len(PcentreL),1):
    h[i] = -1*(np.log(PcentreL[i]/Psurface[i]))*((R*Temp_V1_17[i])/(Mda*g)) # altitude at the centre of level (m)
    h_RH[i] = -1*(np.log(PcentreL[i]/Psurface[i]))*((R*Temp_V1_17[i])/(Mda*g)) # altitude at the centre of level (m)

#------------------------------------------------------------------------------
# APPEND ALTITUDES OF 2M & 10M TO H

# 2d array for 2m
h2m = np.zeros(len(Psurface))
h2m.fill(2)
h2m = np.array([h2m]).T

# 2d array for 10m
h10m = np.zeros(len(Psurface))
h10m.fill(10)
h10m = np.array([h10m]).T

# Append 2m & 10m to H
h = np.concatenate((h10m,h),axis=1)
h = np.concatenate((h2m,h),axis=1)
h = pd.DataFrame.from_dict(h)
h = np.array(h)

#------------------------------------------------------------------------------
# APPEND TEMPERATURES AT 2M & 10M TO TEMP_V1_17

# 2d array for temperature at 2m
T2m = Temp2m_V1_17
T2m = np.array([T2m]).T

# 2d array for temperature at 10m
T10m = Temp10m_V1_17
T10m = np.array([T10m]).T

# Append T2m & T10m to TEMP_V1_17
Temp_V1_17 = np.concatenate((T10m,Temp_V1_17),axis=1)
Temp_V1_17 = np.concatenate((T2m,Temp_V1_17),axis=1)
Temp_V1_17 = pd.DataFrame.from_dict(Temp_V1_17)
Temp_V1_17 = np.array(Temp_V1_17)

#------------------------------------------------------------------------------
# CALCULATE PRESSURE AT 2m & 10m

# variables
Psurface   = SLP_V1_17/100 # the "true" surface pressure (hPa)
g          = 9.80665       # acceleration due to gravity (m/s2)
Mda        = 0.02897       # molar mass dry air (kg/mol)
R          = 8.3143        # specific gas constant (J/K/mol)
T2m        = Temp2m_V1_17  # temperature at 2m (K)
T10m       = Temp10m_V1_17 # temperature at 10m (K)
h2m        = 2             # altitude of 2m
h10m       = 10            # altitude of 10m

# equation
p2m  = Psurface*np.exp(-1*(Mda*g)/(R*Temp2m_V1_17)*h2m)
p10m = Psurface*np.exp(-1*(Mda*g)/(R*Temp10m_V1_17)*h10m)

#------------------------------------------------------------------------------
# APPEND PRESSURE AT 2M & 10M TO PCENTREL

# 2d array for presure at 2m
p2m = np.array([p2m]).T

# 2d array for temperature at 10m
p10m = np.array([p10m]).T

# Append p2m & p10m to PCENTREL
PcentreL = np.concatenate((p10m,PcentreL),axis=1)
PcentreL = np.concatenate((p2m,PcentreL),axis=1)
PcentreL = pd.DataFrame.from_dict(PcentreL)
PcentreL = np.array(PcentreL)

#------------------------------------------------------------------------------
# CALCULATE PRESSURE AT 100m & 1000m

# variables
x = h         # altitude at the centre of level (m)
y = PcentreL  # pressure at the centre of the level (hPa)
y = np.log(y) # need log(pressure) to make linear

# set empty arrays for p100m & p1000m
p100m  = np.zeros([len(Psurface)]) # [Time=1003]
p1000m = np.zeros([len(Psurface)]) # [Time=1003]

# calculate pressure at 100m & 1000m
for i in range(0,len(Psurface),1):
    inter = interpolate.interp1d(x[i],y[i], kind='linear', axis=0)
    #print(inter)
    # need .exp to convert presssure back from log form
    p100m[i]  = np.exp(inter(100))  # pressure at 100m (hPa)
    p1000m[i] = np.exp(inter(1000)) # pressure at 1000m (hPa)

#------------------------------------------------------------------------------
# CALCULATE TEMPERATURE AT 100m & 1000m

# variables
x  = h          # altitude at the centre of level (m)
y2 = Temp_V1_17 # temperature at the centre of the level (K)

# set empty arrays for T100m & T1000m
T100  = np.zeros([len(Psurface)]) # [Time=1003]
T1000 = np.zeros([len(Psurface)]) # [Time=1003]

# calculate temperature at 100m & 1000m
for i in range(0,len(Psurface),1):
    inter2 = interpolate.interp1d(x[i],y2[i], kind='linear', axis=0)
    T100[i]  = inter2(100)  # temperature at 100m (K)
    T1000[i] = inter2(1000) # temperature at 1000m (K)

#------------------------------------------------------------------------------
# Sanity Plot 2: Vertical Interpolation

#I3Test   = xr.open_mfdataset('/Users/ncp532/Documents/Data/MERRA2/V1_17/MERRA2.201*.I3.2x25.nc4', combine='by_coords') # I3 Fields
#TimeTest = I3Test.time[2:1005].values # DateTime (YYYY-MM_DD HH:MM:SS)
#TempTest = I3Test.T[2:1005,:,79,103].values # DateTime (YYYY-MM_DD HH:MM:SS)
#
#fig = plt.figure()
#ax=plt.subplot(111)
#
## Plot Variables
#ax.plot(TempTest[0,:], h[0,:], marker='x', c='blue', markersize = 8.0, linestyle='-', label='Temp Original')
#ax.plot(T100[0], 100, marker='x', c='red', markersize = 8.0, linestyle='-', label='Temp at 100m')
#ax.plot(T1000[0], 1000, marker='x', c='green', markersize = 8.0, linestyle='-', label='Temp at 1000m')
#
## Format x-axis
#ax.tick_params(axis='x',pad=15)
#
## Format y-axis 1
#ax.yaxis.set_major_locator(ticker.MultipleLocator(300))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))
#ax.set_ylim(0,1200)
#
## Plot the axis labels, legend and title
#plt.title('Sanity Test 2: Altitude Interpolation', fontsize=25, y=1.02)
#ax.set_ylabel('Altitude (m)', fontsize=15)
#ax.set_xlabel('Temperature (K)', fontsize=15)
#plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.)

#------------------------------------------------------------------------------
# APPEND EASTWARD WIND AT 2M & 10M TO EWIND_V1_17

# 2d array for eastward wind at 2m
EWind2m = EWind10m_V1_17
EWind2m = np.array([EWind2m]).T

# 2d array for eastward wind at 10m
EWind10m = EWind10m_V1_17
EWind10m = np.array([EWind10m]).T

# Append EWind2m & EWind10m to EWind_V1_17
EWind_V1_17 = np.concatenate((EWind10m,EWind_V1_17),axis=1)
EWind_V1_17 = np.concatenate((EWind2m,EWind_V1_17),axis=1)
EWind_V1_17 = pd.DataFrame.from_dict(EWind_V1_17)
EWind_V1_17 = np.array(EWind_V1_17)

#------------------------------------------------------------------------------
# CALCULATE EASTWARD WIND AT 100m & 1000m

# variables
x  = h           # altitude at the centre of level (m)
y4 = EWind_V1_17 # eastward wind at the centre of the level (m/s)

# set empty arrays for U100m & U1000m
U100m  = np.zeros([len(Psurface)]) # [Time=1003]
U1000m = np.zeros([len(Psurface)]) # [Time=1003]

# calculate eastward wind at 100m & 1000m
for i in range(0,len(Psurface),1):
    inter4 = interpolate.interp1d(x[i],y4[i], kind='linear', axis=0)
    U100m[i]  = inter4(100)  # eastward wind at 100m (m/s)
    U1000m[i] = inter4(1000) # eastward wind at 1000m (m/s)

#------------------------------------------------------------------------------
# APPEND NORTHWARD WIND AT 2M & 10M TO NWIND_V1_17

# 2d array for northward wind at 2m
NWind2m = NWind10m_V1_17
NWind2m = np.array([NWind2m]).T

# 2d array for northward wind at 10m
NWind10m = NWind10m_V1_17
NWind10m = np.array([NWind10m]).T

# Append NWind2m & NWind10m to NWind_V1_17
NWind_V1_17 = np.concatenate((NWind10m,NWind_V1_17),axis=1)
NWind_V1_17 = np.concatenate((NWind2m,NWind_V1_17),axis=1)
NWind_V1_17 = pd.DataFrame.from_dict(NWind_V1_17)
NWind_V1_17 = np.array(NWind_V1_17)
    
#------------------------------------------------------------------------------
# CALCULATE NORTHWARD WIND AT 100m & 1000m

# variables
x  = h           # altitude at the centre of level (m)
y5 = NWind_V1_17 # eastward wind at the centre of the level (m/s)

# set empty arrays for V100m & V1000m
V100m  = np.zeros([len(Psurface)]) # [Time=1003]
V1000m = np.zeros([len(Psurface)]) # [Time=1003]

# calculate northward wind at 100m & 1000m
for i in range(0,len(Psurface),1):
    inter5 = interpolate.interp1d(x[i],y5[i], kind='linear', axis=0)
    V100m[i]  = inter5(100)  # northward wind at 100m (m/s)
    V1000m[i] = inter5(1000) # northward wind at 1000m (m/s)

#------------------------------------------------------------------------------
# CALCULATE WIND SPEED AT 100m & 1000m

WS_V1_17         = np.sqrt(EWind_V1_17**2
                           + NWind_V1_17**2)              # Wind speed (m/s)

WS100m  = np.sqrt(U100m**2 + V100m**2)   # Wind speed at 100m (m/s)
WS1000m = np.sqrt(U1000m**2 + V1000m**2) # Wind speed at 1000m (m/s)

#------------------------------------------------------------------------------
# CALCULATE THE SATURATION VAPOR PRESSURE (es)

# variables
es1    = 6.11          # Reference vapor pressure (hPa)
T1     = 273.15        # Reference temperature (K)
Llv    = 40.8          # Latent heat of vaporisation (kJ/mol) (or 2260 kJ/kg)
R      = 8.3143        # Specific gas constant (J/K/mol)
Mwv    = 0.01802       # Mass water vapor Kg/mol
Rv     = R/Mwv         # Gas constant for water vapor (J/kg/K)
T2m    = Temp2m_V1_17  # Temperature at 10m (K)
T10m   = Temp10m_V1_17 # Temperature at 10m (K)
T100m  = T100          # Temperature at 100m (K)
T1000m = T1000         # Temperature at 1000m (K)
T      = Temp_V1_17    # Temperature at centre of level (K)

# set an empty array for es
es = np.zeros([len(Psurface),74]) # [Time=1003,Levels=72]

# equation
es2m    = es1*np.exp(-1*(Llv/Rv)*((1/T2m)-(1/T1)))    # Saturation vapor pressure at 10m (hPa)
es10m   = es1*np.exp(-1*(Llv/Rv)*((1/T10m)-(1/T1)))   # Saturation vapor pressure at 10m (hPa)
es100m  = es1*np.exp(-1*(Llv/Rv)*((1/T100m)-(1/T1)))  # Saturation vapor pressure at 100m (hPa)
es1000m = es1*np.exp(-1*(Llv/Rv)*((1/T1000m)-(1/T1))) # Saturation vapor pressure at 1000m (hPa)

for i in range(0,len(Psurface),1):
    es[i] = es1*np.exp(-1*(Llv/Rv)*((1/T[i])-(1/T1))) # Saturation vapor pressure at the centre of level (hPa)

#------------------------------------------------------------------------------
# CALCULATE RELATIVE HUMIDITY AT 2m

# variables
QV2m = QV2m_V1_17    # specific humidity at 2m (kg/kg)
p2m  = SLP_V1_17/100 # sea level pressure (hPa)
Mwv  = 0.01802       # molar mass water (kg/mol)
Mda  = 0.02897       # molar mass dry air (kg/mol)
E    = Mwv/Mda       # Ratio of molar mass water to molar mass dry air 
es2m = es2m          # Saturation vapor pressure at 10m (hPa)

# equation
RH2m = (QV2m*p2m)/(E*es2m) 

#------------------------------------------------------------------------------
# CALCULATE RELATIVE HUMIDITY AT 10m

# make array of altitude 2m
h2m = np.zeros(len(Psurface))
h2m.fill(2)
h2m = np.array([h2m]).T

# add the 2m altitude to h_RH
h_RH1 = np.concatenate((h2m,h_RH),axis=1)
h_RH1 = pd.DataFrame.from_dict(h_RH1)
h_RH1 = np.array(h_RH1)

# add the RH at 2m to RH_V1_17
RHT = np.array([RH2m]).T
RHT = np.concatenate((RHT,RH_V1_17),axis=1)
RHT = pd.DataFrame.from_dict(RHT) # relative humidity at 2m (1)
RHT = np.array(RHT)

# variables
x  = h_RH1 # altitude at the centre of level (m)
y3 = RHT   # relative humidity at the centre of the level (1)

# set empty arrays for RH100m & RH1000m
RH10m = np.zeros([len(Psurface)]) # [Time=1003]
RH100m  = np.zeros([len(Psurface)]) # [Time=1003]
RH1000m = np.zeros([len(Psurface)]) # [Time=1003]

# calculate relative humidity at 10m, 100m & 1000m
for i in range(0,len(Psurface),1):
    inter6 = interpolate.interp1d(x[i],y3[i], kind='linear', axis=0)
    RH10m[i]   = inter6(10)   # relative humidity at 10m (1)
    RH100m[i]  = inter6(100)  # relative humidity at 100m (1)
    RH1000m[i] = inter6(1000) # relative humidity at 1000m (1)

#------------------------------------------------------------------------------
# APPEND RELATIVE HUMIDITY AT 2M & 10M TO RH_V1_17

# 2d array for presure at 2m
RH2m = np.array([RH2m]).T

# 2d array for temperature at 10m
RH10m = np.array([RH10m]).T

# Append p2m & p10m to PCENTREL
RH_V1_17 = np.concatenate((RH10m,RH_V1_17),axis=1)
RH_V1_17 = np.concatenate((RH2m,RH_V1_17),axis=1)
RH_V1_17 = pd.DataFrame.from_dict(RH_V1_17)
RH_V1_17 = np.array(RH_V1_17)

#------------------------------------------------------------------------------
# CALCULATE THE WATER VAPOR PRESSURE (e)

# variables
RH2m    = RH2m.flatten()# Relative humidity at 10m (1)
RH10m   = RH10m.flatten()# Relative humidity at 10m (1)
RH100m  = RH100m   # Relative humidity at 100m (1)
RH1000m = RH1000m  # Relative humidity at 1000m (1)
RH      = RH_V1_17 # Relative humidity at centre of level (1)
es10m   = es2m     # Saturation vapor pressure at 2m (hPa)
es10m   = es10m    # Saturation vapor pressure at 10m (hPa)
es100m  = es100m   # Saturation vapor pressure at 100m (hPa)
es1000m = es1000m  # Saturation vapor pressure at 1000m (hPa)
es      = es       # Saturation vapor pressure at the centre of level (hPa)

# set an empty array for e
e = np.zeros([len(Psurface),74]) # [Time=1003,Levels=72]

# equation
e2m    = (RH2m/100)*es2m       # Water vapor pressure at 2m (hPa)
e10m   = (RH10m/100)*es10m     # Water vapor pressure at 10m (hPa)
e100m  = (RH100m/100)*es100m   # Water vapor pressure at 100m (hPa)
e1000m = (RH1000m/100)*es1000m # Water vapor pressure at 1000m (hPa)

for i in range(0,len(Psurface),1):
    e[i] = (RH[i]/100)*es[i] # Water vapor pressure at centre of level (hPa)
    
#------------------------------------------------------------------------------
# CALCULATE THE VIRTUAL TEMPERATURE (Tv)

# variables
T2m    = Temp2m_V1_17  # Temperature at 2m (K)
T10m   = Temp10m_V1_17 # Temperature at 10m (K)
T100m  = T100          # Temperature at 100m (K)
T1000m = T1000         # Temperature at 1000m (K)
T      = Temp_V1_17    # Temperature at centre of level (K)
e2m    = e2m           # Water vapor pressure at 2m (hPa)
e10m   = e10m          # Water vapor pressure at 10m (hPa)
e100m  = e100m         # Water vapor pressure at 100m (hPa)
e1000m = e1000m        # Water vapor pressure at 1000m (hPa)
e      = e             # Water vapor pressure at centre of level (hPa)
p2m    = p2m.flatten() # pressure at 2m (hPa)
p10m   = p10m.flatten()# pressure at 10m (hPa)
p100m  = p100m         # pressure at 100m (hPa)
p1000m = p1000m        # pressure at 1000m (hPa)
p      = PcentreL      # pressure at centre of level (hPa)
Mwv    = 0.01802       # molar mass water (kg/mol)
Mda    = 0.02897       # molar mass dry air (kg/mol)
E      = Mwv/Mda       # Ratio of molar mass water to molar mass dry air 

# set an empty array for Tv
Tv = np.zeros([len(Psurface),74]) # [Time=1003,Levels=72]

# equation
Tv2m    = T2m/(1-((e2m/p2m)*(1-E)))          # Virtual temperature at 2m (K)
Tv10m   = T10m/(1-((e10m/p10m)*(1-E)))       # Virtual temperature at 10m (K)
Tv100m  = T100m/(1-((e100m/p100m)*(1-E)))    # Virtual temperature at 100m (K)
Tv1000m = T1000m/(1-((e1000m/p1000m)*(1-E))) # Virtual temperature at 1000m (K)

for i in range(0,len(Psurface),1):
    Tv[i] = T[i]/(1-((e[i]/p[i])*(1-E))) # Virtual temperature at centre of level (K)
    
#------------------------------------------------------------------------------
# CALCULATE THE VIRTUAL POTENTIAL TEMPERATURE (thetaV)

# variables
Tv2m    = Tv2m     # virtual temperature at 2m (K)
Tv10m   = Tv10m    # virtual temperature at 10m (K)
Tv100m  = Tv100m   # virtual temperature at 100m (K)
Tv1000m = Tv1000m  # virtual temperature at 1000m (K)
p0      = 1000     # reference pressure (hPa) # SHOULD I USE PS_V1_17 or SLP_V1_17 INSTEAD OF 1000?
p2m     = p2m      # pressure at 2m (hPa)
p10m    = p10m     # pressure at 10m (hPa)
p100m   = p100m    # pressure at 100m (hPa)
p1000m  = p1000m   # pressure at 1000m (hPa)
p       = PcentreL # pressure at centre of level (hPa)
R       = 8.314    # specific gas constant (J/K/mol)
Cpd     = 29.07    # specific heat at constant pressure (J/K/mol)

# set an empty array for thetaV
thetaV = np.zeros([len(Psurface),74]) # [Time=1003,Levels=72]

# equation
thetaV2m    = Tv2m*(p0/p2m)**(R/Cpd)       # Virtual potential temperature at 2m (K)
thetaV10m   = Tv10m*(p0/p10m)**(R/Cpd)     # Virtual potential temperature at 10m (K)
thetaV100m  = Tv100m*(p0/p100m)**(R/Cpd)   # Virtual potential temperature at 100m (K)
thetaV1000m = Tv1000m*(p0/p1000m)**(R/Cpd) # Virtual potential temperature at 1000m (K)

for i in range(0,len(Psurface),1):
    thetaV[i] = Tv[i]*(p0/p[i])**(R/Cpd) # Virtual potential temperature at centre of level (K)

#------------------------------------------------------------------------------
# CALCULATE THE AIR DENSITY (rho)
   
# variables
prho = p*100   # pressure at centre of level (Pa = hPa*100)
R    = 8.314   # specific gas constant (J/K/mol)
Mda  = 0.02897 # molar mass dry air (kg/mol)
Rda  = R/Mda   # Gas constant for dry air (J/kg/K)
T    = T       # Temperature at centre of level (K)

# set an empty array for rho
rho = np.zeros([len(Psurface),74]) # [Time=1003,Levels=74]

# equation
for i in range(0,len(Psurface),1):
    rho[i] = prho[i]/(Rda*T[i]) # Air density at centre of level (kg/m3)

# variables
x  = h   # Altitude at the centre of level (m)
y6 = rho # Air density at centre of level (kg/m3)

# set empty arrays for rho100m to rho1900m
rho100m  = np.zeros([len(Psurface)]) # [Time=1003]
rho300m  = np.zeros([len(Psurface)]) # [Time=1003]
rho500m  = np.zeros([len(Psurface)]) # [Time=1003]
rho700m  = np.zeros([len(Psurface)]) # [Time=1003]
rho900m  = np.zeros([len(Psurface)]) # [Time=1003]
rho1100m = np.zeros([len(Psurface)]) # [Time=1003]
rho1300m = np.zeros([len(Psurface)]) # [Time=1003]
rho1500m = np.zeros([len(Psurface)]) # [Time=1003]
rho1700m = np.zeros([len(Psurface)]) # [Time=1003]
rho1900m = np.zeros([len(Psurface)]) # [Time=1003]
rho2100m  = np.zeros([len(Psurface)]) # [Time=1003]
rho2300m  = np.zeros([len(Psurface)]) # [Time=1003]
rho2500m  = np.zeros([len(Psurface)]) # [Time=1003]
rho2700m  = np.zeros([len(Psurface)]) # [Time=1003]
rho2900m  = np.zeros([len(Psurface)]) # [Time=1003]
rho3100m = np.zeros([len(Psurface)]) # [Time=1003]
rho3300m = np.zeros([len(Psurface)]) # [Time=1003]
rho3500m = np.zeros([len(Psurface)]) # [Time=1003]
rho3700m = np.zeros([len(Psurface)]) # [Time=1003]
rho3900m = np.zeros([len(Psurface)]) # [Time=1003]

# calculate air density at 100m to 1900m
for i in range(0,len(Psurface),1):
    inter7      = interpolate.interp1d(x[i],y6[i], kind='linear', axis=0)
    rho100m[i]  = inter7(100)  # Air density at 100m (kg/m3)
    rho300m[i]  = inter7(300)  # Air density at 300m (kg/m3)
    rho500m[i]  = inter7(500)  # Air density at 500m (kg/m3)
    rho700m[i]  = inter7(700)  # Air density at 700m (kg/m3)
    rho900m[i]  = inter7(900)  # Air density at 900m (kg/m3)
    rho1100m[i] = inter7(1100) # Air density at 1100m (kg/m3)
    rho1300m[i] = inter7(1300) # Air density at 1300m (kg/m3)
    rho1500m[i] = inter7(1500) # Air density at 1500m (kg/m3)
    rho1700m[i] = inter7(1700) # Air density at 1700m (kg/m3)
    rho1900m[i] = inter7(1900) # Air density at 1900m (kg/m3)
    rho2100m[i] = inter7(2100) # Air density at 2300m (kg/m3)
    rho2300m[i] = inter7(2300) # Air density at 2300m (kg/m3)
    rho2500m[i] = inter7(2500) # Air density at 2500m (kg/m3)
    rho2700m[i] = inter7(2700) # Air density at 2700m (kg/m3)
    rho2900m[i] = inter7(2900) # Air density at 2900m (kg/m3)
    rho3100m[i] = inter7(3100) # Air density at 3100m (kg/m3)
    rho3300m[i] = inter7(3300) # Air density at 3300m (kg/m3)
    rho3500m[i] = inter7(3500) # Air density at 3500m (kg/m3)
    rho3700m[i] = inter7(3700) # Air density at 3700m (kg/m3)
    rho3900m[i] = inter7(3900) # Air density at 3900m (kg/m3)
    
#------------------------------------------------------------------------------
# CALCULATE THE BULK RICHARDSON NUMBER (Rib)

# variables for calculation
g           = 9.80665     # acceleration due to gravity (m/s2)
z100m       = 100         # altitude of 100m
z1000m      = 1000        # altitude of 1000m
z           = h           # altitude at the centre of level (m)
thetaV10m   = thetaV10m   # virtual potential temperature at surface (K)
thetaV100m  = thetaV100m  # virtual potential temperature at 100m (K)
thetaV1000m = thetaV1000m # virtual potential temperature at 1000m (K)
thetaV      = thetaV      # Virtual potential temperature at centre of level (K)
U100m       = U100m       # Eastward wind at 100m (m/s)
U1000m      = U1000m      # Eastward wind at 1000m (m/s)
U           = EWind_V1_17 # Eastward wind at centre of level (m/s)
V100m       = V100m       # Northward wind at 100m (m/s)
V1000m      = V1000m      # Northward wind at 1000m (m/s)
V           = NWind_V1_17 # Northward wind at centre of level (m/s)

# set an empty array for Rib
Rib = np.zeros([len(Psurface),74]) # [Time=1003,Levels=72]

# equation
Rib_100m  = ((g*z100m)/thetaV10m)*((thetaV100m-thetaV10m)/(U100m**2+V100m**2))     # Bulk Richardson number (100m)
Rib_1000m = ((g*z1000m)/thetaV10m)*((thetaV1000m-thetaV10m)/(U1000m**2+V1000m**2)) # Bulk Richardson number (1000m)

for i in range(0,len(Psurface),1):
    Rib[i] = ((g*z[i])/thetaV10m[i])*((thetaV[i]-thetaV10m[i])/(U[i]**2+V[i]**2)) # Bulk Richardson number at centre of level
    
#------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE PBL HEIGHT

# The PBL height is then chosen as the lowest height at which the value of Rib reaches a critical threshold.
# Swanson et al (2020) uses a critical thresholds of 0.25

def pbl_bulk_richardson( z, thetav, ws, crit=0.25):
# Bulk Richardson definition of PBL height
# z, thetav, and ws are all numpy arrays with values defined vertically 
# (i.e. 2m, 10m, etc.)
    if (max(z) > 100):
        print('pbl_bulk_richardson: z must by in km')
        raise SystemExit()

    g0 = 9.81
    
    # Convert km -> m
    zm = z * 1e3
    # If input is already in meters
    #zm = z
    
    # Bulk Richardson number
    Ri = g0 * zm * ( thetav - thetav[0] ) / thetav[0] / ws**2 
    
    # mixing height is the lowest level where Ri > 0.25
    L=0
    while (Ri[L] < crit ):
        L+=1
    zpbl1 = z[L]
        
    # Use a 2nd order polynomial fit to smooth the profile and 
    # refine the height of the critical value
    pf = np.polyfit( z[L-1:L+2], Ri[L-1:L+2]-crit, 2 )

    # Root from quadratic formula (positive root because Ri is increasing in z)
    zpbl = (-pf[1] + np.sqrt( pf[1]**2 - 4*pf[0]*pf[2] ) ) / (2 * pf[0])
    
    # Add a fraction to L indicating where in the interval the pbl height is
    if (zpbl >= zpbl1):
        L = L + (zpbl-zpbl1) / ( z[L+1] - z[L] )
    else:
        L = L + (zpbl-zpbl1) / ( z[L]   - z[L-1])
    
    return zpbl, L

# variables for function
z      = h/1000   # Altitude at centre of level (km)
thetav = thetaV   # Virtual potential temperature at centre of level (K)
ws     = WS_V1_17 # Wind speed at centre of level (m/s)

# set empty arrays for zpbl & L
zpbl = np.zeros([len(Psurface)]) # height of the PBL (km)
L    = np.zeros([len(Psurface)]) # mixing height (km)

# call the function
for i in range(0,len(Psurface),1):
    zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], crit=0.25)

#------------------------------------------------------------------------------
# BUILD THE NEW DATAFRAMES

# V1_17
dfV1_17 = np.column_stack((PrecTot_V1_17,PrecSno_V1_17,Rainfall_V1_17,SLP_V1_17,PS_V1_17,Temp2m_V1_17,Temp10m_V1_17,T100m,T1000m,EWind10m_V1_17,U100m,U1000m,NWind10m_V1_17,U100m,U1000m,WS10m_V1_17,WS100m,WS1000m,es2m,es10m,es100m,es1000m,e2m,e10m,e100m,e1000m,Tv2m,Tv10m,Tv100m,Tv1000m,thetaV2m,thetaV10m,thetaV100m,thetaV1000m,zpbl,L,
                           rho100m,rho300m,rho500m,rho700m,rho900m,rho1100m,rho1300m,rho1500m,rho1700m,rho1900m,rho2100m,rho2300m,rho2500m,rho2700m,rho2900m,rho3100m,rho3300m,rho3500m,rho3700m,rho3900m))
dfV1_17 = pd.DataFrame.from_dict(dfV1_17)
dfV1_17.columns = ['PrecTot','PrecSno','Rainfall','SLP','SurfPres','Temp2m','Temp10m','Temp100m','Temp1000m','EWind10m','EWind100m','EWind1000m','NWind10m','NWind100m','NWind1000m','WS10m','WS100m','WS1000m','SVP2m','SVP10m','SVP100m','SVP1000m','WVP2m','WVP10m','WVP100m','WVP1000m','VT2m','VT10m','VT100m','VT1000m','VPT2m','VPT10m','VPT100m','VPT1000m','MLH','MLH_FracLev',
                   'rho100m','rho300m','rho500m','rho700m','rho900m','rho1100m','rho1300m','rho1500m','rho1700m','rho1900m','rho2100m','rho2300m','rho2500m','rho2700m','rho2900m','rho3100m','rho3300m','rho3500m','rho3700m','rho3900m']
dfV1_17['DateTime'] = Time_V1_17
dfV1_17 = dfV1_17.set_index('DateTime')

#------------------------------------------------------------------------------
# EXPORT THE DATAFRAMES AS .CSV

dfV1_17.to_csv('/Users/ncp532/Documents/Data/MERRA2/V3_18D_MERRA2.csv')
