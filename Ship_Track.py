#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:19:08 2019

@author: ncp532
"""

# A script that averages observation based on a geos-chem grid
# The same way GEOS-Chem averages the planeflight output
# This code is only valid for the 2x25 grid! Not the 4x5 grid!
# Test it with the plane.log files if you don't trust me :D 
import numpy as np
import csv
from datetime import datetime, timedelta
import time
from functools import reduce
#-------------------------------------------------------------------------
#Theoretically you need to change 3 things 
#1. TimeStep 
#2. Year of the measurement 
#3. Path to your data

#However you might additionally need to modify:
#4. The way you read in the data based on your measurement file
#5. Additionally eliminate some weird values from file if needed
#6. Specify the directory where to write the new file

# 1. Specify the timestep
TimeStep = 20
# 2. Year of the measurements, to convert DOY back to normal time
year = 2012
# 3. Specify your path
path = '/home/bb907/Desktop/Computer/Work/Data/Ship_Train_insitu/Ship_2012_2013/Measurements/'

#The name of your measurement file
data=path+'Measurements/SouthernSurveyor'+str(year)+'_Ooofti.dat'

#Later add this:
#Specify the grid (type 2.5 for 2x2.5 or 5 for 4x5)

#-------------------------------------------------------------------------
#Information about the grid 
#this is for : CTM_TYPE('GEOS5_47L', RES=2)
nlon       = 144
nlat       = 91
lonmid     = np.arange(-180, 180, 2.5)
lonedge    = np.arange(-181.250, 180, 2.5)
latmid     = np.arange(-90.0, 91.0, 2) 
#because of the poles the first and last array will be -89.5 and 89.5 and not -90 and 90
latmid[0]  =-89.5
latmid[-1] = 89.5
latedge    =np.arange(-91.0, 92.0, 2)
#same for the edges, it will be -90 and 90, and not -91, 91
latedge[0] =-90
latedge[-1]= 90
#Because this is for the ship data, the levels and pressures are not important here
#lev=1, it is all on the surface

#-------------------------------------------------------------------------

#These are the variables that we need from the file
date = []       #date of measurements
lat  = []       #latitude of measurements
lon  = []       #longitude of measurements
pres = []       #pressure of measurements
co2  = []       #the co2 measurement
ch4  = []       #the ch4 measurement  
co   = []       #the co measurement

#4. Read in the data
with open(data) as csvfile:
    read=csv.reader(csvfile, delimiter=' ')
    for i in range(31):    #to skip the header
        tmp=read.next()
      #  print tmpex
    for row in read:
        date.append(row[0])
        lat.append(float(row[1]))
        lon.append(float(row[2]))
        pres.append(float(row[3]))
        co2.append(float(row[4]))
        ch4.append(float(row[6]))
        co.append(float(row[7]))

#5. Eliminate the wrong values (the 99999.99)
goodlat=np.where([(s != 99999.99) and (s!= 99999.99999) for s in lat])[0]
goodlon=np.where([(s != 99999.99) and (s!= 99999.99999) for s in lon])[0]
goodpres=np.where([(s != 99999.99) and (s!= 99999.99999) and (s!= 99999.9) for s in pres])[0]
goodco2=np.where([(s != 99999.99) and (s!= 99999.99999) for s in co2])[0]
goodch4=np.where([(s != 99999.99) and (s!= 99999.99999) for s in ch4])[0]
goodco=np.where([(s != 99999.99) and (s!= 99999.99999) for s in co])[0]

#Get the good indicies where I have all the values
#Use reduce if I want the intersect of more than 2 arrays
goodvals = reduce(np.intersect1d, (goodlat,goodlon,goodpres, goodco2, goodch4, goodco))
#for all the data take only the lines with all good values
date=[date[i] for i in goodvals]
lat1=["%6.2f"%(lat[i]) for i in goodvals]   #I need this formating so that it is the same as the GC output
lon1=["%6.2f"%(lon[i]) for i in goodvals]   #Otherwise if the decimal values are different it will mess up the averaging
pres=[pres[i] for i in goodvals]
co2=[co2[i] for i in goodvals]
ch4=[ch4[i] for i in goodvals]
co=[co[i] for i in goodvals]

#make them float
lat =[]
lon = []
for i in range(len(lat1)):
    lat.append(float(lat1[i]))
    lon.append(float(lon1[i]))

#The measurements are seconds from some date, convert it to normal time       
time0=datetime(year,1,1,0)
dates= [time0+timedelta(seconds=int(dd)) for dd in date]

#Convert time to string
datestr = [datetime.strftime(d,"%d/%m/%Y %H:%M:%S") for d in dates]

#Calculate Day of the year, use Julian day? with decimal numbers
#This should be OK, it gives the same numbers as the IDL scripts, wuhuuu
DOY = np.zeros(len(datestr))
for i in range(len(datestr)):
        modh = dates[i].hour
        modm = dates[i].minute   
        #print date[i]
        DOY[i] = time.strptime(datestr[i], "%d/%m/%Y %H:%M:%S").tm_yday+modh/24.+modm/(60.0*24)
        #DOY = zenaz.zenaz(sitelat,sitelon,sitealt,modyy,modmm,moddd,modfday)

# Find the index of the lat, lon, and level for each observation
ctm_Lat_INDEX1 = ( np.interp( lat, latmid, range(nlat)) )
ctm_Lon_INDEX1 = ( np.interp( lon, lonmid, range(nlon) ) )
#I HAD A PROBLEM HERE, because if the number is 24.5 IDL rounds it as 25 while python rounds it as 24 so it was a mismtach...
#so don't use numpy.round use the round function, so round it up
ctm_Lat_INDEX=[]
ctm_Lon_INDEX=[]
for i in range(len(ctm_Lat_INDEX1)):
    ctm_Lat_INDEX.append(round(ctm_Lat_INDEX1[i]))
    ctm_Lon_INDEX.append(round(ctm_Lon_INDEX1[i]))

#create the arrays 
ctm_Lat_INDEX=np.array(ctm_Lat_INDEX)
ctm_Lon_INDEX=np.array(ctm_Lon_INDEX)
 
#In IDL interpol extrapolates beyond the ends of the indgen array so we need to cap the INDEX
#But in python np.interp is not doing that, so we are OK, no need to fix that

# Now assign a unique integer to each time interval of width minutes
# Time_INDEX_full gives the number of TimeStep intervals that have
# elapsed since the start of the year
Time_INDEX_full = np.int32(DOY * 24. * 60. / TimeStep )

#We can reduce the magnitude of Time_INDEX_full, by starting
#Time_INDEX at 0 for the first observation
Time_INDEX = Time_INDEX_full - min( Time_INDEX_full )
nTime = max( Time_INDEX )+1

# Now construct a single number to describe the 4D grouping
Group_INDEX = np.int32( Time_INDEX * nlat * nlon ) +( ctm_Lat_INDEX * nlon ) +( ctm_Lon_INDEX )

#Create a function that will average different arrays
def tapply(array, group, function):
    #Find the unique values in the classification array GROUP
    groupvalues=np.unique(group)
    # Make an array with the same type as the input array
    result = np.zeros( len( groupvalues ))
    # Loop over the number of unique values
    for i in range(len(groupvalues)):
       # Find which elements share the same value
       index = np.where( group == groupvalues[i])[0]
       # make numpy arrays from the input arrays
       array = np.array(array)
       # Apply the given function to the common elements
       result[i] = function(array[index])
    return result

#The averaged values
Lat_geos = tapply(lat, Group_INDEX, np.nanmean)
Lon_geos = tapply( lon, Group_INDEX, np.nanmean )
Prs_geos = tapply( pres, Group_INDEX, np.nanmean)
DOY_geos = tapply( DOY, Group_INDEX, np.nanmean)
CO2_geos = tapply( co2, Group_INDEX, np.nanmean)
CH4_geos = tapply( ch4, Group_INDEX, np.nanmean)
CO_geos = tapply( co, Group_INDEX, np.nanmean)

#Convert DOY back to normal time
time0=datetime(year,1,1,0)
date_geos= [time0+timedelta(days=d) for d in DOY_geos]

#datetime counts from 1, so there is an offset of 1 day when I convert the DOY
#If I add the DOY, it will add the days to 1 (and not 0)
delta = timedelta(days= 1) #offset
offset = []
for i in date_geos:
    offset.append(i-delta)
date_geos = offset

#Separate the date into date and time
Date_geos = []
Time_geos = []
for i in range(len(date_geos)):
    Date_geos.append(date_geos[i].strftime("%Y%m%d"))
    Time_geos.append(date_geos[i].strftime("%H%M"))

#-------------------------------------------------------------------------
#6. Write the result into a new file
with open('avg_meas'+str(TimeStep)+'min'+str(year)+'shipAll.csv', 'wb') as f:
   writer = csv.DictWriter(f, delimiter = ' ', fieldnames = ['Date', 'Time', 'DOY','Lat','Lon', 'Press', 'CO2', 'CH4', 'CO'])
   writer.writeheader()
   for i in range(len(DOY_geos)):
       writer.writerow({'Date':Date_geos[i], 'Time':Time_geos[i], 'DOY':"%6.3f"%(DOY_geos[i]), 'Lat':"%6.3f"%(Lat_geos[i]),
       'Lon':"%6.3f"%(Lon_geos[i]),'Press':"%6.2f"%(Prs_geos[i]), 'CO2':"%6.2f"%(CO2_geos[i]), 'CH4':"%7.2f"%(CH4_geos[i]),
       'CO':"%5.2f"%(CO_geos[i])})
       
#The below part is just a specific way of writing things into files...       
#Based on the date create seperate file for every day
       
#Get only the days when we have measurements, and sort them (only days)
days =list( set([x for x in Date_geos]))
#Sort them because set messes up the order
days.sort(key=lambda x: time.mktime(time.strptime(x,"%Y%m%d")))

# Separate the file into smaller files by date
for day in days:
    #dayindexes=np.where(datum == day)[0]`
    dayindexes= np.where([ d == day for d in Date_geos ])[0]
    dayy= [ Date_geos[d] for d in dayindexes ]
    Year=str(day[:4])
    Month=str(day[4:6])
    Day=str(day[6:])
    with open(path+'Avg_Measurements/'+Year+'_sorted/avg_meas'+str(TimeStep)+'min'+Year+'_'+Month+'_'+Day+'_ship.csv', 'wb') as f:
       writer = csv.DictWriter(f, delimiter = ' ', fieldnames = ['Date', 'Time', 'DOY','Lat','Lon', 'Press', 'CO2', 'CH4', 'CO'])
       writer.writeheader()
       for i in dayindexes:
           writer.writerow({'Date':Date_geos[i], 'Time':Time_geos[i], 'DOY':"%6.3f"%(DOY_geos[i]), 'Lat':"%6.3f"%(Lat_geos[i]),
           'Lon':"%6.3f"%(Lon_geos[i]),'Press':"%6.2f"%(Prs_geos[i]), 'CO2':"%6.2f"%(CO2_geos[i]), 'CH4':"%7.2f"%(CH4_geos[i]),
           'CO':"%5.2f"%(CO_geos[i])})
       
