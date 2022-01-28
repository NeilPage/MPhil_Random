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

# Radiosonde
#RS_2017   = Dataset('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/Radiosonde/ARM_Radiosonde_2017.nc')
RS_2017 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/Radiosonde/ARM_Radiosonde_20171101_20180324.csv',   index_col=0)

# Met
#Met_V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ShipTrack/V1_17_underway_60.csv', index_col=0)
#Met_V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ShipTrack/V2_17_underway_60.csv', index_col=0)
#Met_V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ShipTrack/V3_17_underway_60.csv', index_col=0)

#------------------------------------------------------------------------------
# SET THE DATE

# Radiosonde
RS_2017.index   = (pd.to_datetime(RS_2017.index,   dayfirst=True))

# Met
#Met_V1_17.index = (pd.to_datetime(Met_V1_17.index, dayfirst=True))
#Met_V2_17.index = (pd.to_datetime(Met_V2_17.index, dayfirst=True))
#Met_V3_17.index = (pd.to_datetime(Met_V3_17.index, dayfirst=True))

#------------------------------------------------------------------------------
# RESAMPLE THE MET DATASETS TO 1 HOUR TIME RESOLUTION
#
#Met_V1_17 = Met_V1_17.resample('60T').mean()
#Met_V2_17 = Met_V2_17.resample('60T').mean()
#Met_V3_17 = Met_V3_17.resample('60T').mean()
#
##------------------------------------------------------------------------------
## COMBINE THE MET DATASETS
#
#Met_2017  = pd.concat([Met_V1_17,Met_V2_17,Met_V3_17], axis =0)

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

## Radiosonde
#DateTime  = RS_2017['DateTime']  # DateTime ()
#Altitude  = pd.DataFrame(np.array(RS_2017['altitude']))  # Altitude (m)
#Pressure  = RS_2017.['pressure'][:,:]  # Pressure (hPa)
#TempDry   = RS_2017['tempdry']   # Temperature (C)
#WindSpeed = RS_2017['windspeed'] # Wind speed (m/s)
#WindDir   = RS_2017['winddir']   # Wind direction (deg)
#RelHum    = RS_2017['relhum']    # Relative humidity (%)
#EastWind  = RS_2017['eastwind']  # Eastward wind component (m/s)
#NorthWind = RS_2017['northwind'] # Northward wind component (m/s)

#------------------------------
# Variables
#------------------------------
Time  = RS_2017.index        # UT Time (YYYY-MM-DD HH:MM:SS)
Alt   = RS_2017['alt']       # Altitude (m)
Pres  = RS_2017['pres']      # Pressure (hPa)
Tdry  = RS_2017['tdry']      # Dry bulb temperature (C)
Tdry  = Tdry + 273.15        # Dry bulb temperature (K)
Wspd  = RS_2017['wspd']      # Wind speed (m/s)
Wdir  = RS_2017['deg']       # Wind direction (deg)
RH    = RS_2017['rh']        # Relative humidity (%)
Uwind = RS_2017['u_wind']    # Eastward wind component (m/s)
Vwind = RS_2017['v_wind']    # Northward wind component (m/s)

#------------------------------
# Seperate the individual sondes
#------------------------------
 # array of sub-arrays (starts with first value)
UTTime    = [[Time[0]]]
Altitude  = [[Alt[0]]]
Pressure  = [[Pres[0]]]
TempDry   = [[Tdry[0]]]
WindSpeed = [[Wspd[0]]]
WindDir   = [[Wdir[0]]]
RelHum    = [[RH[0]]]
EastWind  = [[Uwind[0]]]
NorthWind = [[Vwind[0]]]

# go through each element based on the altitude
for i in range(1, len(Alt)):
    # If altitude is larger than previous altitude
    if Alt[i - 1] <= Alt[i]:
        # Add the variables to the last sub-array
        Altitude[len(Altitude)   - 1].append(Alt[i])
        Pressure[len(Pressure)   - 1].append(Pres[i])
        UTTime[len(UTTime)       - 1].append(Time[i])
        TempDry[len(TempDry)     - 1].append(Tdry[i])
        WindSpeed[len(WindSpeed) - 1].append(Wspd[i])
        WindDir[len(WindDir)     - 1].append(Wdir[i])
        RelHum[len(RelHum)       - 1].append(RH[i])
        EastWind[len(EastWind)   - 1].append(Uwind[i])
        NorthWind[len(NorthWind) - 1].append(Vwind[i])
    # If the altitude isn't larger than previous altitude
    else:
        # Add the variables to a new sub-array
        Altitude.append([Alt[i]])
        Pressure.append([Pres[i]])
        UTTime.append([Time[i]])
        TempDry.append([Tdry[i]])
        WindSpeed.append([Wspd[i]])
        WindDir.append([Wdir[i]])
        RelHum.append([RH[i]])
        EastWind.append([Uwind[i]])
        NorthWind.append([Vwind[i]])

# Convert the results to Pandas DataFrames 
Altitude  = pd.DataFrame(Altitude)
Pressure  = pd.DataFrame(Pressure)
UTTime    = pd.DataFrame(UTTime)
TempDry   = pd.DataFrame(TempDry)
WindSpeed = pd.DataFrame(WindSpeed)
WindDir   = pd.DataFrame(WindDir)
RelHum    = pd.DataFrame(RelHum)
EastWind  = pd.DataFrame(EastWind)
NorthWind = pd.DataFrame(NorthWind)

# Select a subset of the Dataframes
Altitude  = (Altitude.T[0:400]).T
Pressure  = (Pressure.T[0:400]).T
UTTime    = (UTTime.T[0:400]).T
TempDry   = (TempDry.T[0:400]).T
WindSpeed = (WindSpeed.T[0:400]).T
WindDir   = (WindDir.T[0:400]).T
RelHum    = (RelHum.T[0:400]).T
EastWind  = (EastWind.T[0:400]).T
NorthWind = (NorthWind.T[0:400]).T
 
# Convert UTTime to an index for the other DataFrames
IndexTime = UTTime.T.loc[0]

# Levels
Levels = np.array(UTTime.columns)

# Set the DataFrame index to IndexTime
UTTime.index    = IndexTime
Altitude.index  = IndexTime
Pressure.index  = IndexTime
TempDry.index   = IndexTime
WindSpeed.index = IndexTime
WindDir.index   = IndexTime
RelHum.index    = IndexTime
EastWind.index  = IndexTime
NorthWind.index = IndexTime

# Resample the dataset to 1-hour resolution
Altitude  = Altitude.resample('60T').mean()
Pressure  = Pressure.resample('60T').mean()
TempDry   = TempDry.resample('60T').mean()
WindSpeed = WindSpeed.resample('60T').mean()
WindDir   = WindDir.resample('60T').mean()
RelHum    = RelHum.resample('60T').mean()
EastWind  = EastWind.resample('60T').mean()
NorthWind = NorthWind.resample('60T').mean()

# Interpolate the values within the dataset
Altitude  = Altitude.interpolate(method='time')
Pressure  = Pressure.interpolate(method='time')
TempDry   = TempDry.interpolate(method='time')
WindSpeed = WindSpeed.interpolate(method='time')
WindDir   = WindDir.interpolate(method='time')
RelHum    = RelHum.interpolate(method='time')
EastWind  = EastWind.interpolate(method='time')
NorthWind = NorthWind.interpolate(method='time')
IndexTime = Altitude.index

#------------------------------------------------------------------------------
# CALCULATE PRESSURE AT 100m, 200m & 1000m

# variables
x = np.array(Altitude.reset_index(drop=True)) # altitude at the centre of level (m)
y = np.array(Pressure .reset_index(drop=True)) # # pressure at the centre of the level (hPa)
y = np.log(y) # need log(pressure) to make linear

# set empty arrays for p100m & p1000m
p100m  = np.zeros([len(IndexTime)]) # [Time=1003]
p200m  = np.zeros([len(IndexTime)]) # [Time=1003]
p1000m = np.zeros([len(IndexTime)]) # [Time=1003]

# calculate pressure at 100m & 1000m
for i in range(0,len(IndexTime),1):
    inter = interpolate.interp1d(x[i],y[i], kind='linear', axis=0)
    #print(inter)
    # need .exp to convert presssure back from log form
    p100m[i]  = np.exp(inter(100))  # pressure at 100m (hPa)
    p200m[i]  = np.exp(inter(200))  # pressure at 200m (hPa)
    p1000m[i] = np.exp(inter(1000)) # pressure at 1000m (hPa)

#------------------------------------------------------------------------------
# CALCULATE TEMPERATURE AT 100m, 200m & 1000m

# variables
x  = np.array(Altitude.reset_index(drop=True)) # altitude at the centre of level (m)
y2 = np.array(TempDry.reset_index(drop=True))  # temperature at the centre of the level (K)

# set empty arrays for T100m & T1000m
T100  = np.zeros([len(IndexTime)]) # [Time=1003]
T200  = np.zeros([len(IndexTime)]) # [Time=1003]
T1000 = np.zeros([len(IndexTime)]) # [Time=1003]

# calculate temperature at 100m & 1000m
for i in range(0,len(IndexTime),1):
    inter2 = interpolate.interp1d(x[i],y2[i], kind='linear', axis=0)
    T100[i]  = inter2(100)  # temperature at 100m (K)
    T200[i]  = inter2(200)  # temperature at 200m (K)
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
# CALCULATE EASTWARD WIND AT 100m & 1000m

# variables
x  = np.array(Altitude.reset_index(drop=True)) # altitude at the centre of level (m)
y4 = np.array(EastWind.reset_index(drop=True)) # eastward wind at the centre of the level (m/s)

# set empty arrays for U100m & U1000m
U100m  = np.zeros([len(IndexTime)]) # [Time=1003]
U200m  = np.zeros([len(IndexTime)]) # [Time=1003]
U1000m = np.zeros([len(IndexTime)]) # [Time=1003]

# calculate eastward wind at 100m & 1000m
for i in range(0,len(IndexTime),1):
    inter4 = interpolate.interp1d(x[i],y4[i], kind='linear', axis=0)
    U100m[i]  = inter4(100)  # eastward wind at 100m (m/s)
    U200m[i]  = inter4(200)  # eastward wind at 100m (m/s)
    U1000m[i] = inter4(1000) # eastward wind at 1000m (m/s)

#------------------------------------------------------------------------------
# CALCULATE NORTHWARD WIND AT 100m & 1000m

# variables
x  = np.array(Altitude.reset_index(drop=True)) # altitude at the centre of level (m)
y5 = np.array(NorthWind.reset_index(drop=True)) # northward wind at the centre of the level (m/s)

# set empty arrays for V100m & V1000m
V100m  = np.zeros([len(IndexTime)]) # [Time=1003]
V200m  = np.zeros([len(IndexTime)]) # [Time=1003]
V1000m = np.zeros([len(IndexTime)]) # [Time=1003]

# calculate northward wind at 100m & 1000m
for i in range(0,len(IndexTime),1):
    inter5 = interpolate.interp1d(x[i],y5[i], kind='linear', axis=0)
    V100m[i]  = inter5(100)  # northward wind at 100m (m/s)
    V200m[i]  = inter5(200)  # northward wind at 100m (m/s)
    V1000m[i] = inter5(1000) # northward wind at 1000m (m/s)

#------------------------------------------------------------------------------
# CALCULATE WIND SPEED AT 100m & 1000m

WS = np.sqrt(EastWind**2 + NorthWind**2)              # Wind speed (m/s)

WS100m  = np.sqrt(U100m**2  + V100m**2)  # Wind speed at 100m (m/s)
WS200m  = np.sqrt(U200m**2  + V200m**2)  # Wind speed at 200m (m/s)
WS1000m = np.sqrt(U1000m**2 + V1000m**2) # Wind speed at 1000m (m/s)

#------------------------------------------------------------------------------
# CALCULATE THE SATURATION VAPOR PRESSURE (es)

# variables
es1    = 6.11              # Reference vapor pressure (hPa)
T1     = 273.15            # Reference temperature (K)
Llv    = 40.8              # Latent heat of vaporisation (kJ/mol) (or 2260 kJ/kg)
R      = 8.3143            # Specific gas constant (J/K/mol)
Mwv    = 0.01802           # Mass water vapor Kg/mol
Rv     = R/Mwv             # Gas constant for water vapor (J/kg/K)
T17m   = np.array(TempDry) # Temperature at 17m (K)
T17m   = T17m[:,0]         # Temperature at 17m (K)
T100m  = T100              # Temperature at 100m (K)
T200m  = T200              # Temperature at 200m (K)
T1000m = T1000             # Temperature at 1000m (K)
T      = np.array(TempDry.reset_index(drop=True)) # Temperature at centre of level (K)

# set an empty array for es
es = np.zeros([len(IndexTime),400]) # [Time=1003,Levels=72]

# equation
es17m   = es1*np.exp(-1*(Llv/Rv)*((1/T17m)-(1/T1)))   # Saturation vapor pressure at 10m (hPa)
es100m  = es1*np.exp(-1*(Llv/Rv)*((1/T100m)-(1/T1)))  # Saturation vapor pressure at 100m (hPa)
es200m  = es1*np.exp(-1*(Llv/Rv)*((1/T200m)-(1/T1)))  # Saturation vapor pressure at 10m (hPa)
es1000m = es1*np.exp(-1*(Llv/Rv)*((1/T1000m)-(1/T1))) # Saturation vapor pressure at 1000m (hPa)

for i in range(0,len(IndexTime),1):
    es[i] = es1*np.exp(-1*(Llv/Rv)*((1/T[i])-(1/T1))) # Saturation vapor pressure at the centre of level (hPa)

#--------c----------------------------------------------------------------------
# CALCULATE RELATIVE HUMIDITY

# variables
x  = np.array(Altitude.reset_index(drop=True))   # altitude at the centre of level (m)
y3 = np.array(RelHum.reset_index(drop=True))/100 # relative humidity at the centre of the level (1)
RH17m   = np.array(RelHum)/100                   # relative humidity at 17m (1)
RH17m   = RH17m[:,0]                             # relative humidity at 17m (1)

# set empty arrays for RH100m & RH1000m
RH100m  = np.zeros([len(IndexTime)]) # [Time=1003]
RH200m  = np.zeros([len(IndexTime)]) # [Time=1003]
RH1000m = np.zeros([len(IndexTime)]) # [Time=1003]

# calculate relative humidity at 10m, 100m & 1000m
for i in range(0,len(IndexTime),1):
    inter6 = interpolate.interp1d(x[i],y3[i], kind='linear', axis=0)
    RH100m[i]  = inter6(100)  # relative humidity at 100m (1)
    RH200m[i]  = inter6(200)  # relative humidity at 200m (1)
    RH1000m[i] = inter6(1000) # relative humidity at 1000m (1)

#------------------------------------------------------------------------------
# CALCULATE THE WATER VAPOR PRESSURE (e)

# variables
RH17m   = RH17m    # Relative humidity at 17m (1)
RH100m  = RH100m   # Relative humidity at 100m (1)
RH200m  = RH200m   # Relative humidity at 200m (1)
RH1000m = RH1000m  # Relative humidity at 1000m (1)
RH      = np.array(RelHum.reset_index(drop=True)) # relative humidity at the centre of the level (1)
es17m   = es17m    # Saturation vapor pressure at 17m (hPa)
es100m  = es100m   # Saturation vapor pressure at 100m (hPa)
es200m  = es200m   # Saturation vapor pressure at 200m (hPa)
es1000m = es1000m  # Saturation vapor pressure at 1000m (hPa)
es      = es       # Saturation vapor pressure at the centre of level (hPa)

# set an empty array for e
e = np.zeros([len(IndexTime),400]) # [Time=1003,Levels=72]

# equation
e17m   = (RH17m/100)*es17m     # Water vapor pressure at 17m (hPa)
e100m  = (RH100m/100)*es100m   # Water vapor pressure at 100m (hPa)
e200m  = (RH200m/100)*es200m   # Water vapor pressure at 200m (hPa)
e1000m = (RH1000m/100)*es1000m # Water vapor pressure at 1000m (hPa)

for i in range(0,len(IndexTime),1):
    e[i] = (RH[i]/100)*es[i] # Water vapor pressure at centre of level (hPa)
    
#------------------------------------------------------------------------------
# CALCULATE THE VIRTUAL TEMPERATURE (Tv)

# variables
T17m   = T17m               # Temperature at 17m (K)
T100m  = T100m              # Temperature at 100m (K)
T200m  = T200m              # Temperature at 200m (K)
T1000m = T1000m             # Temperature at 1000m (K)
T      = T                  # Temperature at centre of level (K)
e17m   = e17m               # Water vapor pressure at 17m (hPa)
e100m  = e100m              # Water vapor pressure at 100m (hPa)
e200m  = e200m              # Water vapor pressure at 200m (hPa)
e1000m = e1000m             # Water vapor pressure at 1000m (hPa)
e      = e                  # Water vapor pressure at centre of level (hPa)
p17m   = np.array(Pressure) # pressure at 17m (hPa)
p17m   = p17m[:,0]          # presure at 17m (hPa)
p100m  = p100m              # pressure at 100m (hPa)
p200m  = p200m              # pressure at 200m (hPa)
p1000m = p1000m             # pressure at 1000m (hPa)
p      = np.array(Pressure .reset_index(drop=True)) # pressure at centre of level (hPa)
Mwv    = 0.01802            # molar mass water (kg/mol)
Mda    = 0.02897            # molar mass dry air (kg/mol)
E      = Mwv/Mda            # Ratio of molar mass water to molar mass dry air 

# set an empty array for Tv
Tv = np.zeros([len(IndexTime),400]) # [Time=1003,Levels=72]

# equation
Tv17m   = T17m/(1-((e17m/p17m)*(1-E)))       # Virtual temperature at 17m (K)
Tv100m  = T100m/(1-((e100m/p100m)*(1-E)))    # Virtual temperature at 100m (K)
Tv200m  = T200m/(1-((e200m/p200m)*(1-E)))    # Virtual temperature at 200m (K)
Tv1000m = T1000m/(1-((e1000m/p1000m)*(1-E))) # Virtual temperature at 1000m (K)

for i in range(0,len(IndexTime),1):
    Tv[i] = T[i]/(1-((e[i]/p[i])*(1-E))) # Virtual temperature at centre of level (K)
    
#------------------------------------------------------------------------------
# CALCULATE THE VIRTUAL POTENTIAL TEMPERATURE (thetaV)

# variables
Tv17m   = Tv17m    # virtual temperature at 17m (K)
Tv100m  = Tv100m   # virtual temperature at 100m (K)
Tv200m  = Tv200m   # virtual temperature at 200m (K)
Tv1000m = Tv1000m  # virtual temperature at 1000m (K)
p0      = 1000     # reference pressure (hPa) # SHOULD I USE PS_V1_17 or SLP_V1_17 INSTEAD OF 1000?
p17m    = p17m     # pressure at 17m (hPa)
p100m   = p100m    # pressure at 100m (hPa)
p200m   = p200m    # pressure at 200m (hPa)
p1000m  = p1000m   # pressure at 1000m (hPa)
p       = p        # pressure at centre of level (hPa)
R       = 8.314    # specific gas constant (J/K/mol)
Cpd     = 29.07    # specific heat at constant pressure (J/K/mol)

# set an empty array for thetaV
thetaV = np.zeros([len(IndexTime),400]) # [Time=1003,Levels=72]

# equation
thetaV17m   = Tv17m*(p0/p17m)**(R/Cpd)     # Virtual potential temperature at 17m (K)
thetaV100m  = Tv100m*(p0/p100m)**(R/Cpd)   # Virtual potential temperature at 100m (K)
thetaV200m  = Tv200m*(p0/p200m)**(R/Cpd)   # Virtual potential temperature at 200m (K)
thetaV1000m = Tv1000m*(p0/p1000m)**(R/Cpd) # Virtual potential temperature at 1000m (K)

for i in range(0,len(IndexTime),1):
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
rho = np.zeros([len(IndexTime),400]) # [Time=1003,Levels=74]

# equation
for i in range(0,len(IndexTime),1):
    rho[i] = prho[i]/(Rda*T[i]) # Air density at centre of level (kg/m3)

# variables
x  = np.array(Altitude.reset_index(drop=True)) # altitude at the centre of level (m)
y6 = rho # Air density at centre of level (kg/m3)

# set empty arrays for rho100m to rho1900m
rho100m  = np.zeros([len(IndexTime)]) # [Time=1003]
rho300m  = np.zeros([len(IndexTime)]) # [Time=1003]
rho500m  = np.zeros([len(IndexTime)]) # [Time=1003]
rho700m  = np.zeros([len(IndexTime)]) # [Time=1003]
rho900m  = np.zeros([len(IndexTime)]) # [Time=1003]
rho1100m = np.zeros([len(IndexTime)]) # [Time=1003]
rho1300m = np.zeros([len(IndexTime)]) # [Time=1003]
rho1500m = np.zeros([len(IndexTime)]) # [Time=1003]

# calculate air density at 100m to 1900m
for i in range(0,len(IndexTime),1):
    inter7      = interpolate.interp1d(x[i],y6[i], kind='linear', axis=0)
    rho100m[i]  = inter7(100)  # Air density at 100m (kg/m3)
    rho300m[i]  = inter7(300)  # Air density at 300m (kg/m3)
    rho500m[i]  = inter7(500)  # Air density at 500m (kg/m3)
    rho700m[i]  = inter7(700)  # Air density at 700m (kg/m3)
    rho900m[i]  = inter7(900)  # Air density at 900m (kg/m3)
    rho1100m[i] = inter7(1100) # Air density at 1100m (kg/m3)
    # rho1300m[i] = inter7(1300) # Air density at 1300m (kg/m3)
    # rho1500m[i] = inter7(1500) # Air density at 1500m (kg/m3)
    
#------------------------------------------------------------------------------
# CALCULATE THE BULK RICHARDSON NUMBER (Rib)

# variables for calculation
g           = 9.80665     # acceleration due to gravity (m/s2)
z100m       = 100         # altitude of 100m
z200m       = 200         # altitude of 200m
z1000m      = 1000        # altitude of 1000m
z           = np.array(Altitude.reset_index(drop=True)) # altitude at the centre of level (m)
thetaV17m   = thetaV17m   # virtual potential temperature at surface (K)
thetaV100m  = thetaV100m  # virtual potential temperature at 100m (K)
thetaV200m  = thetaV200m  # virtual potential temperature at 200m (K)
thetaV1000m = thetaV1000m # virtual potential temperature at 1000m (K)
thetaV      = thetaV      # Virtual potential temperature at centre of level (K)
U100m       = U100m       # Eastward wind at 100m (m/s)
U200m       = U100m       # Eastward wind at 200m (m/s)
U1000m      = U1000m      # Eastward wind at 1000m (m/s)
U           = np.array(EastWind.reset_index(drop=True))   # Eastward wind at centre of level (m/s)
U17m        = U[:,0]      # Eastward wind at 17m (m/s)
V100m       = V100m       # Northward wind at 100m (m/s)
V200m       = V200m       # Northward wind at 200m (m/s)
V1000m      = V1000m      # Northward wind at 1000m (m/s)
V           = np.array(NorthWind.reset_index(drop=True)) # Northward wind at centre of level (m/s)
V17m        = V[:,0]      # Northward wind at 17m (m/s)
WS17m       = np.sqrt(U17m**2 + V17m**2)  # Wind speed at 17m (m/s)

# set an empty array for Rib
Rib = np.zeros([len(IndexTime),400]) # [Time=1003,Levels=72]

# equation
Rib_100m  = ((g*z100m)/thetaV17m)*((thetaV100m-thetaV17m)/(U100m**2+V100m**2))     # Bulk Richardson number (100m)
Rib_200m  = ((g*z200m)/thetaV17m)*((thetaV200m-thetaV17m)/(U200m**2+V200m**2))     # Bulk Richardson number (200m)
Rib_1000m = ((g*z1000m)/thetaV17m)*((thetaV1000m-thetaV17m)/(U1000m**2+V1000m**2)) # Bulk Richardson number (1000m)

for i in range(0,len(IndexTime),1):
    Rib[i] = ((g*z[i])/thetaV17m[i])*((thetaV[i]-thetaV17m[i])/(U[i]**2+V[i]**2)) # Bulk Richardson number at centre of level
    
#------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE PBL HEIGHT

# The PBL height is then chosen as the lowest height at which the value of Rib reaches a critical threshold.
# Swanson et al (2020) uses a critical thresholds of 0.25

#def pbl_bulk_richardson( z, thetav, ws, Rib, crit=0.25):
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
    #Ri = Rib
    Ri = g0 * zm * ( thetav - thetav[0] ) / thetav[0] / ws**2 
    #print(Ri[0])
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
z      = np.array(Altitude.reset_index(drop=True))/1000   # Altitude at centre of level (km)
thetav = thetaV   # Virtual potential temperature at centre of level (K)
#ws     = np.sqrt(U**2 + V**2) # wind speed (m/s)
ws     = np.array(WindSpeed.reset_index(drop=True)) # Wind speed at centre of level (m/s)

# pbl 1
Time1   = IndexTime[0:1206]
z1      = z[0:1206]
thetav1 = thetaV[0:1206]
ws1     = ws[0:1206]

# pbl 2
Time2   = IndexTime[1207:3439]
z2      = z[1207:3439]
thetav2 = thetaV[1207:3439]
ws2     = ws[1207:3439]

# set empty arrays for zpbl & L
zpbl = np.zeros([len(IndexTime)]) # height of the PBL (km)
L    = np.zeros([len(IndexTime)]) # mixing height (km)
Ri   = np.zeros([len(IndexTime)]) # mixing height (km)

# zpbl1
zpbl1 = np.zeros([len(z1)]) # height of the PBL (km)
L1    = np.zeros([len(z1)]) # mixing height (km)
Ri1   = np.zeros([len(z1)]) # mixing height (km)

# zpbl1
zpbl2 = np.zeros([len(z2)]) # height of the PBL (km)
L2    = np.zeros([len(z2)]) # mixing height (km)
Ri2   = np.zeros([len(z2)]) # mixing height (km)

# call the function
# for i in range(0,len(IndexTime),1):
# #    zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], Rib[i], crit=0.25)
#     zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], crit=0.25)

for i in range(0,len(z1),1):
#    zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], Rib[i], crit=0.25)
#    zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], crit=0.25)
    zpbl1[i], L1[i] = pbl_bulk_richardson(z1[i], thetav1[i], ws1[i], crit=0.25)

for i in range(0,len(z2),1):
#    zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], Rib[i], crit=0.25)
#    zpbl[i], L[i] = pbl_bulk_richardson(z[i], thetav[i], ws[i], crit=0.25)
    zpbl2[i], L2[i] = pbl_bulk_richardson(z2[i], thetav2[i], ws2[i], crit=0.25)

#------------------------------------------------------------------------------
# Concat zpbl1 & zpbl2

zpblNaN = np.nan

test = np.append(zpbl1,[zpblNaN])

# Append zpblNaN to zpbl1
zpbl = np.concatenate((test,zpbl2))

#------------------------------------------------------------------------------
# BUILD THE NEW DATAFRAMES

# V1_17
dfV1_17 = np.column_stack((p17m,p100m,p200m,p1000m,T17m,T100m,T200m,T1000m,U17m,U100m,U200m,U1000m,V17m,V100m,V200m,V1000m,WS17m,WS100m,WS200m,WS1000m,es17m,es100m,es200m,es1000m,e17m,e100m,e200m,e1000m,Tv17m,Tv100m,Tv200m,Tv1000m,thetaV17m,thetaV100m,thetaV200m,thetaV1000m,zpbl,L,
                           rho100m,rho300m,rho500m,rho700m,rho900m,rho1100m))
dfV1_17 = pd.DataFrame.from_dict(dfV1_17)
dfV1_17.columns = ['Pres17m','Pres100m','Pres200m','Pres1000m','Temp17m','Temp100m','Temp200m','Temp1000m','EWind17m','EWind100m','EWind200m','EWind1000m','NWind17m','NWind100m','NWind200m','NWind1000m','WS17m','WS100m','WS200m','WS1000m','SVP17m','SVP100m','SVP200m','SVP1000m','WVP17m','WVP100m','WVP200m','WVP1000m','VT17m','VT100m','VT200m','VT1000m','VPT17m','VPT100m','VPT200m','VPT1000m','MLH','MLH_FracLev',
                   'rho100m','rho300m','rho500m','rho700m','rho900m','rho1100m']
dfV1_17['DateTime'] = IndexTime
dfV1_17 = dfV1_17.set_index('DateTime')
dfV1_17 = dfV1_17.dropna()

#------------------------------------------------------------------------------
# EXPORT THE DATAFRAMES AS .CSV

dfV1_17.to_csv('/Users/ncp532/Documents/Data/MERRA2/CAMMPCAN_2017_Radiosonde_Interp1hr.csv')
