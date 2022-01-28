 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 16:21:15 2019

@author: ncp532
"""

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
from windrose import WindroseAxes

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# CAMMPCAN 2017-18
V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V1_17_Data.csv',header=0,encoding = 'unicode_escape')
V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V2_17_Data.csv',header=0,encoding = 'unicode_escape')
V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V3_17_Data.csv',header=0,encoding = 'unicode_escape')

# CAMMPCAN 2018-19
V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V1_18_Data.csv',header=0,encoding = 'unicode_escape')
V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V2_18_Data.csv',header=0,encoding = 'unicode_escape')
V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V3_18_Data.csv',header=0,encoding = 'unicode_escape')

# SIPEXII 2012
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/SIPEXII_Data.csv',header=0,encoding = 'unicode_escape')

#------------------------------------------------------------------------------
# Set the date

V1_17['DateTime'] = pd.to_datetime(V1_17['DateTime']) # Davis timezone is UT+7
V2_17['DateTime'] = pd.to_datetime(V2_17['DateTime']) # Casey timezone is UT+8
V3_17['DateTime'] = pd.to_datetime(V3_17['DateTime']) # Mawson timezone is UT+5
V1_18['DateTime'] = pd.to_datetime(V1_18['DateTime']) # Davis timezone is UT+7
V2_18['DateTime'] = pd.to_datetime(V2_18['DateTime']) # Casey timezone is UT+8
V3_18['DateTime'] = pd.to_datetime(V3_18['DateTime']) # Mawson timezone is UT+5
SIPEXII['DateTime'] = pd.to_datetime(SIPEXII['DateTime']) # SIPEXII timezone is UT+5

#------------------------------------------------------------------------------
# Filter the datasets for midday hours only

#-----------------------------
# CAMMPCAN 2017-18
#-----------------------------
# V1_17 Davis (07:00 to 18:00)
start_time = '07:00:00'
end_time = '18:00:00'
Midday = (V1_17['Time'] >= start_time) & (V1_17['Time'] < end_time)
V1_17_MM = V1_17[Midday]

# V2_17 Casey (08:00 to 16:00)
start_time = '08:00:00'
end_time = '16:00:00'
Midday = (V2_17['Time'] >= start_time) & (V2_17['Time'] < end_time)
V2_17_MM = V2_17[Midday]

# V3_17 Mawson (08:00 to 18:00)
start_time = '08:00:00'
end_time = '18:00:00'
Midday = (V3_17['Time'] >= start_time) & (V3_17['Time'] < end_time)
V3_17_MM = V3_17[Midday]

#-----------------------------
# CAMMPCAN 2018-19
#-----------------------------
# V1_18 Davis (07:00 to 18:00)
start_time = '07:00:00'
end_time = '18:00:00'
Midday = (V1_18['Time'] >= start_time) & (V1_18['Time'] < end_time)
V1_18_MM = V1_18[Midday]

# V2_18 Casey (08:00 to 16:00)
start_time = '08:00:00'
end_time = '16:00:00'
Midday = (V2_18['Time'] >= start_time) & (V2_18['Time'] < end_time)
V2_18_MM = V2_18[Midday]

# V3_18 Mawson (08:00 to 18:00)
start_time = '08:00:00'
end_time = '18:00:00'
Midday = (V3_18['Time'] >= start_time) & (V3_18['Time'] < end_time)
V3_18_MM = V3_18[Midday]

#-----------------------------
# SIPEXII 2012
#-----------------------------
# SIPEXII (07:00 to 18:00)
start_time = '07:00:00'
end_time = '18:00:00'
Midday = (SIPEXII['Time'] >= start_time) & (SIPEXII['Time'] < end_time)
SIPEXII_MM = SIPEXII[Midday]

#------------------------------------------------------------------------------
# Filter dataframe for when filter is less than 60%

V1_17F = (V1_17_MM['Filter'] < 0.6)
V1_17T = V1_17_MM[V1_17F]

V2_17F = (V2_17_MM['Filter'] < 0.6)
V2_17T = V2_17_MM[V2_17F]

V3_17F = (V3_17_MM['Filter'] < 0.6)
V3_17T = V3_17_MM[V3_17F]

V1_18F = (V1_18_MM['Filter'] < 0.6)
V1_18T = V1_18_MM[V1_18F]

V2_18F = (V2_18_MM['Filter'] < 0.6)
V2_18T = V2_18_MM[V2_18F]

V3_18F = (V3_18_MM['Filter'] < 0.6)
V3_18T = V3_18_MM[V3_18F]

SIPEXIIF = (SIPEXII_MM['Filter'] < 0.6)
SIPEXIIT = SIPEXII_MM[SIPEXIIF]

#------------------------------------------------------------------------------
# Define the variables


#------------------------------------------------------------------------------
# Define the variables

# CAMMPCAN (2017-18)
BrO_V1_17 = np.array(V1_17T['surf_vmr(ppmv)']) * 1e6
Lat_V1_17 = np.array(V1_17T['LATITUDE'])
Long_V1_17 = np.array(V1_17T['LONGITUDE'])
WD_s1_17 = np.array(V1_17T['WND_DIR_STRBD_CORR_DEG'])
WD_p1_17 = np.array(V1_17T['WND_DIR_PORT_CORR_DEG'])
WS_s1_17 = np.array(V1_17T['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_p1_17 = np.array(V1_17T['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_V1_17 = (WS_s1_17 + WS_p1_17)/2
WD_vectV1_17 = ((WD_s1_17 * WS_s1_17) / (WS_s1_17 + WS_p1_17)) + ((WD_p1_17 * WS_p1_17) / (WS_s1_17 + WS_p1_17))

BrO_V2_17 = np.array(V2_17T['surf_vmr(ppmv)']) * 1e6
Lat_V2_17 = np.array(V2_17T['LATITUDE'])
Long_V2_17 = np.array(V2_17T['LONGITUDE'])
WD_s2_17 = np.array(V2_17T['WND_DIR_STRBD_CORR_DEG'])
WD_p2_17 = np.array(V2_17T['WND_DIR_PORT_CORR_DEG'])
WS_s2_17 = np.array(V2_17T['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_p2_17 = np.array(V2_17T['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_V2_17 = (WS_s2_17 + WS_p2_17)/2
WD_vectV2_17 = ((WD_s2_17 * WS_s2_17) / (WS_s2_17 + WS_p2_17)) + ((WD_p2_17 * WS_p2_17) / (WS_s2_17 + WS_p2_17))

BrO_V3_17 = np.array(V3_17T['surf_vmr(ppmv)']) * 1e6
Lat_V3_17 = np.array(V3_17T['LATITUDE'])
Long_V3_17 = np.array(V3_17T['LONGITUDE'])
WD_s3_17 = np.array(V3_17T['WND_DIR_STRBD_CORR_DEG'])
WD_p3_17 = np.array(V3_17T['WND_DIR_PORT_CORR_DEG'])
WS_s3_17 = np.array(V3_17T['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_p3_17 = np.array(V3_17T['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_V3_17 = (WS_s3_17 + WS_p3_17)/2
WD_vectV3_17 = ((WD_s3_17 * WS_s3_17) / (WS_s3_17 + WS_p3_17)) + ((WD_p3_17 * WS_p3_17) / (WS_s3_17 + WS_p3_17))

# CAMMPCAN (2018-19)
BrO_V1_18 = np.array(V1_18T['surf_vmr(ppmv)']) * 1e6
Lat_V1_18 = np.array(V1_18T['LATITUDE'])
Long_V1_18 = np.array(V1_18T['LONGITUDE'])
WD_s1_18 = np.array(V1_18T['WND_DIR_STRBD_CORR_DEG'])
WD_p1_18 = np.array(V1_18T['WND_DIR_PORT_CORR_DEG'])
WS_s1_18 = np.array(V1_18T['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_p1_18 = np.array(V1_18T['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_V1_18 = (WS_s1_18 + WS_p1_18)/2
WD_vectV1_18 = ((WD_s1_18 * WS_s1_18) / (WS_s1_18 + WS_p1_18)) + ((WD_p1_18 * WS_p1_18) / (WS_s1_18 + WS_p1_18))

BrO_V2_18 = np.array(V2_18T['surf_vmr(ppmv)']) * 1e6
Lat_V2_18 = np.array(V2_18T['LATITUDE'])
Long_V2_18 = np.array(V2_18T['LONGITUDE'])
WD_s2_18 = np.array(V2_18T['WND_DIR_STRBD_CORR_DEG'])
WD_p2_18 = np.array(V2_18T['WND_DIR_PORT_CORR_DEG'])
WS_s2_18 = np.array(V2_18T['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_p2_18 = np.array(V2_18T['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_V2_18 = (WS_s2_18 + WS_p2_18)/2
WD_vectV2_18 = ((WD_s2_18 * WS_s2_18) / (WS_s2_18 + WS_p2_18)) + ((WD_p2_18 * WS_p2_18) / (WS_s2_18 + WS_p2_18))

BrO_V3_18 = np.array(V3_18T['surf_vmr(ppmv)']) * 1e6
Lat_V3_18 = np.array(V3_18T['LATITUDE'])
Long_V3_18 = np.array(V3_18T['LONGITUDE'])
WD_s3_18 = np.array(V3_18T['WND_DIR_STRBD_CORR_DEG'])
WD_p3_18 = np.array(V3_18T['WND_DIR_PORT_CORR_DEG'])
WS_s3_18 = np.array(V3_18T['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_p3_18 = np.array(V3_18T['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_V3_18 = (WS_s3_18 + WS_p3_18)/2
WD_vectV3_18 = ((WD_s3_18 * WS_s3_18) / (WS_s3_18 + WS_p3_18)) + ((WD_p3_18 * WS_p3_18) / (WS_s3_18 + WS_p3_18))

# SIPEXII (2012)
BrO_SIPEXII = np.array(SIPEXIIT['surf_vmr(ppmv)']) * 1e6
Lat_SIPEXII = np.array(SIPEXIIT['LATITUDE'])
Long_SIPEXII = np.array(SIPEXIIT['LONGITUDE'])
WD_sSIPEXII = np.array(SIPEXIIT['WND_DIR_STRBD_CORR_DEG'])
WD_pSIPEXII = np.array(SIPEXIIT['WND_DIR_PORT_CORR_DEG'])
WS_sSIPEXII = np.array(SIPEXIIT['WND_SPD_STRBD_CORR_KNOT']) * 0.514444444
WS_pSIPEXII = np.array(SIPEXIIT['WND_SPD_PORT_CORR_KNOT']) * 0.514444444
WS_SIPEXII = (WS_sSIPEXII + WS_pSIPEXII)/2
WD_vectSIPEXII = ((WD_sSIPEXII * WS_sSIPEXII) / (WS_sSIPEXII + WS_pSIPEXII)) + ((WD_pSIPEXII * WS_pSIPEXII) / (WS_sSIPEXII + WS_pSIPEXII))

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# SET THE DATE AND TIME
#------------------------------------
# V1_17
dat1 = np.array(V1_17T['Date'])
tim1 = np.array(V1_17T['Time'])
dattim1 = dat1+' '+tim1

#CONVERT TO DATETIME FROM STRING
date1=[]
for i in range(len(dattim1)):
    date1.append(datetime.strptime(dattim1[i],'%d/%m/%Y %H:%M:%S')) # midday data    

#------------------------------------    
# V2_17
dat2 = np.array(V2_17T['Date'])
tim2 = np.array(V2_17T['Time'])
dattim2 = dat2+' '+tim2

#CONVERT TO DATETIME FROM STRING
date2=[]
for i in range(len(dattim2)):
    date2.append(datetime.strptime(dattim2[i],'%d/%m/%Y %H:%M:%S')) # midday data 

#------------------------------------
# V3_17
dat3 = np.array(V3_17T['Date'])
tim3 = np.array(V3_17T['Time'])
dattim3 = dat3+' '+tim3

#CONVERT TO DATETIME FROM STRING
date3=[]
for i in range(len(dattim3)):
    date3.append(datetime.strptime(dattim3[i],'%d/%m/%Y %H:%M:%S')) # midday data 

#------------------------------------
# V1_18
dat4 = np.array(V1_18T['Date'])
tim4 = np.array(V1_18T['Time'])
dattim4 = dat4+' '+tim4

#CONVERT TO DATETIME FROM STRING
date4=[]
for i in range(len(dattim4)):
    date4.append(datetime.strptime(dattim4[i],'%d/%m/%Y %H:%M:%S')) # midday data 

#------------------------------------
# V2_18
dat5 = np.array(V2_18T['Date'])
tim5 = np.array(V2_18T['Time'])
dattim5 = dat5+' '+tim5

#CONVERT TO DATETIME FROM STRING
date5=[]
for i in range(len(dattim5)):
    date5.append(datetime.strptime(dattim5[i],'%d/%m/%Y %H:%M:%S')) # midday data 

#------------------------------------
# V3_18
dat6 = np.array(V3_18T['Date'])
tim6 = np.array(V3_18T['Time'])
dattim6 = dat6+' '+tim6

#CONVERT TO DATETIME FROM STRING
date6=[]
for i in range(len(dattim6)):
    date6.append(datetime.strptime(dattim6[i],'%d/%m/%Y %H:%M:%S')) # midday data 

#------------------------------------
# SIPEXII
dat7 = np.array(SIPEXIIT['Date'])
tim7 = np.array(SIPEXIIT['Time'])
dattim7 = dat7+' '+tim7

#CONVERT TO DATETIME FROM STRING
date7=[]
for i in range(len(dattim7)):
    date7.append(datetime.strptime(dattim7[i],'%d/%m/%Y %H:%M:%S')) # midday data 

#------------------------------------------------------------------------------
# BUILD A NEW DATAFRAME FOR IMPORT TO R AND SAVE AS A .CSV 

dfV1_17 = np.column_stack((dat1, tim1, WD_vectV1_17, WS_V1_17, BrO_V1_17))
dfV1_17 = pd.DataFrame.from_dict(dfV1_17)
dfV1_17.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfV1_17.to_csv('/Users/ncp532/Documents/Data/V1_17.csv')

dfV2_17 = np.column_stack((dat2, tim2, WD_vectV2_17, WS_V2_17, BrO_V2_17))
dfV2_17 = pd.DataFrame.from_dict(dfV2_17)
dfV2_17.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfV2_17.to_csv('/Users/ncp532/Documents/Data/V2_17.csv')

dfV3_17 = np.column_stack((dat3, tim3, WD_vectV3_17, WS_V3_17, BrO_V3_17))
dfV3_17 = pd.DataFrame.from_dict(dfV3_17)
dfV3_17.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfV3_17.to_csv('/Users/ncp532/Documents/Data/V3_17.csv')

dfV1_18 = np.column_stack((dat4, tim4, WD_vectV1_18, WS_V1_18, BrO_V1_18))
dfV1_18 = pd.DataFrame.from_dict(dfV1_18)
dfV1_18.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfV1_18.to_csv('/Users/ncp532/Documents/Data/V1_18.csv')

dfV2_18 = np.column_stack((dat5, tim5, WD_vectV2_18, WS_V2_18, BrO_V2_18))
dfV2_18 = pd.DataFrame.from_dict(dfV2_18)
dfV2_18.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfV2_18.to_csv('/Users/ncp532/Documents/Data/V2_18.csv')

dfV3_18 = np.column_stack((dat6, tim6, WD_vectV3_18, WS_V3_18, BrO_V3_18))
dfV3_18 = pd.DataFrame.from_dict(dfV3_18)
dfV3_18.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfV3_18.to_csv('/Users/ncp532/Documents/Data/V3_18.csv')

dfSIPEXII = np.column_stack((dat7, tim7, WD_vectSIPEXII, WS_SIPEXII, BrO_SIPEXII))
dfSIPEXII = pd.DataFrame.from_dict(dfSIPEXII)
dfSIPEXII.columns = ['Date','Time','WD_vector','Wind_Speed','BrO_(pptv)']
dfSIPEXII.to_csv('/Users/ncp532/Documents/Data/SIPEXII.csv')