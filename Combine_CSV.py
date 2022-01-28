#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:38:52 2019

@author: ncp532
"""

# Data handing packages
import os
import glob
import pandas as pd
import numpy as np

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time
#------------------------------------------------------------------------------

# Set the location for the working directory
os.chdir("/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/Processed_CO_N2O/")

# Set the file type to .csv
extension = 'csv'
# Save a list of all the file names to the variable all_filenames 
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

# Combine all files in the list
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])

#------------------------------------------------------------------------------
# SET THE DATE AND TIME

dattim = np.array(combined_csv['Date'])

#CONVERT TO DATETIME FROM STRING
date=[]
for i in range(len(dattim)):
    date.append(datetime.strptime(dattim[i],'%Y-%m-%d %H:%M:%S'))

combined_csv['date'] = date
combined_csv.set_index('date',inplace=True)
O3 = combined_csv.resample('1S').mean()
##------------------------------------------------------------------------------
## Define the variables
#CO = np.array(combined_csv['CO']) # CO
#N2O = np.array(combined_csv['N2O']) # N2O
#CO_Dry = np.array(combined_csv['CO_Dry']) # CO_Dry
#N2O_Dry = np.array(combined_csv['N2O_Dry']) # N2O_Dry
#
## Calculate the 10 min mean
#def onemin(x, date):
#    df = pd.DataFrame({'X':x}, index=date) 
#    df = df.resample('1T').mean()
#    #Reset the index
#    df =df.reset_index()
#    #extract the values
#    x=df['X']
#    date=df['index']  
#    #convert the pandas series date to list
#    date = date.tolist()
#    return x,date 
## Resample from 1 second to 1 minute resolution   
##datim_1min, date_1min=onemin(dattim,date) # dattim
#CO_1min, date_1min=onemin(CO,date) # CO
#N2O_1min, date_1min=onemin(N2O,date) # N2O
#CO_Dry_1min, date_1min=onemin(CO_Dry,date) # CO_Dry
#N2O_Dry_1min, date_1min=onemin(N2O_Dry,date) # N2O_Dry
##------------------------------------------------------------------------------
## BUILD A NEW DATAFRAME FOR IMPORT TO R AND SAVE AS A .CSV 
#dfARM = np.column_stack((date_1min,CO_1min,N2O_1min,CO_Dry_1min,N2O_Dry_1min))
#dfARM = pd.DataFrame.from_dict(dfARM)
#dfARM.columns = ['Date','CO','N2O','CO_Dry','N2O_Dry']

# Export the combined .csv
O3.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/ARM/all_CO_1sec.csv')
#dfARM.to_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/dfARM.csv')
