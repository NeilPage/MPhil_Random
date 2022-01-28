#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:21:42 2020

@author: ncp532
"""

# DATA HANDLING PACKAGES
import numpy as np
import pandas as pd
import os

#-------------------------------
# RETRIEVE ICE CONTACT TIME

# Set the working directory
os.chdir("/Users/ncp532/Documents/data/SeaIce_Trajectories/2000m")

# Set the folder location
data_folder = '/Users/ncp532/Documents/data/SeaIce_Trajectories/2000m/'

# Function to get all files within the folder
def get_all_files(directory, extension='.csv'):
    dir_list = os.listdir(directory)
    csv_files = []
    for e in dir_list:
        if e.endswith(extension):
            csv_files.append(os.path.realpath(e))
    return csv_files

# Save the files to a list
csv_files = get_all_files(data_folder)

# Set empty arrays for the variables
Year     = np.zeros(len(csv_files))
Month    = np.zeros(len(csv_files))
Day      = np.zeros(len(csv_files))
Hour     = np.zeros(len(csv_files))
Min      = np.zeros(len(csv_files))
Ice_100m = np.zeros(len(csv_files))
Ice_MLH  = np.zeros(len(csv_files))
SIC_100m = np.zeros(len(csv_files))
SIC_MLH = np.zeros(len(csv_files))

# Retrieve the variables
for i in range(len(csv_files)):
    test        = pd.read_csv(csv_files[i])
    Year[i]     = test['Traj Year'][0]
    Month[i]    = test['Traj Mon'][0]
    Day[i]      = test['Traj Day'][0]
    Hour[i]     = test['Traj Hour'][0]
    Min[i]      = test['Traj Min'][0]
    Ice_100m[i] = np.sum(test['Traj over Sea Ice and height < 100 m?'])
    Ice_MLH[i]  = np.sum(test['Traj over Sea Ice and height < mixed layer top?'])
    SIC_100m[i] = np.sum(test['Traj over Sea Ice and height < 100 m?']*test['Sea Ice Conc (0-1)'])
    SIC_MLH[i] = np.sum(test['Traj over Sea Ice and height < mixed layer top?']*test['Sea Ice Conc (0-1)'])
    
#------------------------------------------------------------------------------
# BUILD THE NEW DATAFRAMES

# IceContact
dfIceContact = np.column_stack((Year,Month,Day,Hour,Min,Ice_100m,Ice_MLH,SIC_100m,SIC_MLH))
dfIceContact = pd.DataFrame.from_dict(dfIceContact)
dfIceContact.columns = ['Year','Month','Day','Hour','Min','Ice_100m','Ice_MLH','SIC_100m','SIC_MLH']

# Set the datetime
dfIceContact['DateTime'] = pd.to_datetime(dfIceContact[["Year", "Month", "Day", "Hour"]])
dfIceContact.sort_values('DateTime')
dfIceContact = dfIceContact.set_index('DateTime')

#------------------------------------------------------------------------------
# EXPORT THE DATAFRAMES AS .CSV

dfIceContact.to_csv('/Users/ncp532/Documents/Data/SeaIce_Trajectories/IceContactTime_2000m.csv')
    