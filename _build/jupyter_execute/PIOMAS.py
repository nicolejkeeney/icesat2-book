# Historical ice thickness

This notebook plots ICESat-2 sea ice thickness vs. PIOMAS sea ice thickness. 
 
**Input**:
 - NETCDF4 file from the jupyter book's google bucket
 
 
 **Output**: 
  - Plots of sea ice thickness 

```{tip}
Try running this notebook in Google Colab! Toggle over the rocketship icon at the top of the page and click Colab to open a new window and run the notebook. <br><br>To run a single cell, type **Shift+Enter**. To run the whole notebook, under **Runtime** click **Run all**. Note that you will have to run the notebook from the very beginning and load all the Google Colab dependencies for the code to work.
```

#this cell will load dependencies for running the notebook in Google Colab
#this cell may take a while to run
import sys

#if code is running in google colab, run these cells to install neccessary libraries
if 'google.colab' in sys.modules: 
    !apt-get install -qq libgdal-dev libproj-dev
    !pip install --no-binary shapely shapely --force
    !pip install -q pyproj
    !pip install cartopy
    !pip install netcdf4
    !pip install xarray==0.16.0

## Import notebook dependencies

import os
import numpy as np
import xarray as xr
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt

# Ignore warnings in the notebook to improve display
# You might want to remove this when debugging/writing new code
import warnings
warnings.filterwarnings('ignore')
 
#increase resolution for notebook outputs
%matplotlib inline 
%config InlineBackend.figure_format = 'retina' 

#import utils function 
if 'utils.py' not in os.listdir(os.getcwd()): 
    !gsutil cp gs://icesat2-book-data/utils.py ./
import utils

## Load data into notebook
Copy file from the book's google bucket and load into an xarray dataset. 

### Load ICESat-2 dataset

!gsutil -m cp gs://icesat2-book-data/icesat2-book-winter-data.nc ./
IS2dataset = xr.open_dataset('icesat2-book-winter-data.nc')

### Load PIOMAS dataset

!gsutil -m cp gs://icesat2-book-data/piomas-regridded-data.nc ./
PIOdataset = xr.open_dataset('piomas-regridded-data.nc')

## Restrict datasets to region of interest

#myKeys = [10,11,12,13,15] #Inner Arctic
myKeys = [15] 
IS2dataset = utils.restrictRegionally(IS2dataset, regionKeyList = myKeys)
PIOdataset = utils.restrictRegionally(PIOdataset, regionKeyList = myKeys)

## ICESat-2 mean monthly ice thickness & uncertainty

#calculate mean monthly ICESat-2 sea ice thickness 
IS2means = IS2dataset['ice_thickness_filled'].mean(dim = ['x','y'], skipna = True)
IS2means.attrs['long_name'] = 'ICESat-2 mean ice thickness'

#calculate mean monthly ICESat-2 sea ice thickness uncertainty
IS2uncs = IS2dataset['ice_thickness_unc_filled'].mean(dim = ['x','y'], skipna = True)
IS2uncs.attrs['long_name'] = 'ICESat-2 mean uncertainty'

#### Compile into one dataset

IS2 = xr.Dataset(data_vars = {'mean_ice_thickness': IS2means, 'mean_ice_thickness_unc': IS2uncs})
IS2.attrs = IS2dataset.attrs

## PIOMAS mean monthly ice thickness & uncertainty

### Mean monthly ice thickness

PIOmeans = PIOdataset['PIOMAS_ice_thickness'].mean(dim = ['x','y'], skipna = True)
PIOmeans.attrs['long_name'] = 'PIOMAS mean ice thickness'

### Mean detrended uncertainty

#### Define functions used in calculation 

def get_varDT(Years, Thickness):
    """ Detrend linear time series  
    """
    trendT, interceptT, r_valsT, probT, stderrT = scipy.stats.linregress(Years, Thickness)
    lineT = (trendT * Years) + interceptT
    ThicknessDT = Thickness - lineT
    return ThicknessDT, lineT

def is_mon(xarray_val, desired_mon):
    """ Check if current month is desired month """
    xarray_mon = pd.to_datetime(xarray_val).month
    if xarray_mon == desired_mon:
        return xarray_val
    else: 
        return None

#### Call functions to perform uncertainty calculation

#loop through months in a year and calculate uncertainty for each month 
standDevDT, standDev, monDates, thicknessDT, slopes = [],[],[],[],[]
months = np.arange(1, 12 + 1, 1)  #integers corresponding to months of the year

for mon in months:
    #get dates for each month (i.e. for January get 1978-01, 1979-01, 1980-01, etc)
    times_by_mon = [is_mon(val, mon) for val in PIOmeans.time.values if (is_mon(val, mon) != None)]
    monDates.append(times_by_mon)
    
    #get array with year-mon converted to integer value
    time_integers = np.arange(0, len(times_by_mon))
    
    #get detrended uncertainty and append to list
    thickness_detrended, regLine = get_varDT(time_integers, PIOmeans.sel(time = times_by_mon).values)
    thicknessDT.append(thickness_detrended)
    standDevDT.append(np.std(thickness_detrended))
    slopes.append(regLine)
    
    #get trended uncertainty and append to list 
    standDev.append(np.std(PIOmeans.sel(time = times_by_mon).values))
    
#get a list as long as PIOmeans of all the monthly natural variabilities
standDevCyclic = []
for date in PIOmeans.time.values:
    standDevCyclic.append(standDevDT[pd.to_datetime(date).month - 1])
    
#create DataArray with corresponding uncertainties
PIOuncs = xr.DataArray(standDevCyclic, dims = 'time', coords = {'time': PIOmeans.time.values}, attrs = {'long_name':'PIOMAS mean detrended uncertainty'})

#### Compile into one dataset 

PIO = xr.Dataset(data_vars = {'mean_ice_thickness': PIOmeans, 'mean_ice_thickness_unc': PIOuncs})
PIO.attrs = PIOdataset.attrs

## Winter comparison: PIOMAS vs. ICESat-2
Here we will plot PIOMAS and ICESat-2 mean monthly ice thickness and uncertainty on the same axes to visualize differences in the two datasets over winter 2018-2019 and winter 2019-2020. 

### Define plotting function 

def plotMeansWithUnc(PIO, IS2, title = None):  
    
    #initialize figure and axes
    fig = plt.figure(figsize = [8,5]) 
    ax = plt.axes()
    gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) #add gridlines 

    #plot PIOMAS data 
    PIO['mean_ice_thickness'].plot(ax = ax, color = 'purple', marker = 'o', label = PIO['mean_ice_thickness'].attrs['long_name'])
    ax.fill_between(PIO.time.values, PIO['mean_ice_thickness'] + PIO['mean_ice_thickness_unc'], PIO['mean_ice_thickness'] - PIO['mean_ice_thickness_unc'], facecolor = 'purple', alpha = 0.1, edgecolor = 'none', label = PIO['mean_ice_thickness_unc'].attrs['long_name'])

    #plot ICESat-2 data
    IS2['mean_ice_thickness'].plot(color = 'dodgerblue', marker = 's', linestyle = '-.', label = IS2['mean_ice_thickness'].attrs['long_name'])
    ax.fill_between(IS2.time.values, IS2['mean_ice_thickness'] + IS2['mean_ice_thickness_unc'], IS2['mean_ice_thickness'] - IS2['mean_ice_thickness_unc'], facecolor = 'dodgerblue', alpha = 0.1, edgecolor = 'none', label = IS2['mean_ice_thickness_unc'].attrs['long_name'])

    #display plot with descriptive legend, ylabel, and title
    plt.legend()
    plt.ylabel('Mean sea ice thickness (m)')
    plt.title(PIO.attrs['regions with data'] + ' mean monthly sea ice thickness')
    plt.show()

### Plot comparison for winter 18-19

winter_18_19 = utils.getWinterDateRange(2018,2019)
plotMeansWithUnc(PIO = PIO.sel(time = winter_18_19), IS2 = IS2.sel(time = winter_18_19))

winter_19_20 = utils.getWinterDateRange(2019,2020)
plotMeansWithUnc(PIO = PIO.sel(time = winter_19_20), IS2 = IS2.sel(time = winter_19_20))

## Plot PIOMAS with overlayed ICESat-2 data
From the PIOMAS reanalysis data, this plot shows the yearly melting and growing of sea ice as well as the long term decreasing trend in mean ice thickness. ICESat-2 data is overlayed on the PIOMAS plot, along with a rolling 12 month mean. 

#initialize figure and axes 
fig = plt.figure(figsize = [8,5]) 
ax = plt.axes()
gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) #add gridlines 

#plot data 
PIO['mean_ice_thickness'].plot(ax = ax, color = 'steelblue', label = PIO['mean_ice_thickness'].attrs['long_name'])
IS2['mean_ice_thickness'].plot(ax = ax, color = 'magenta', label = IS2['mean_ice_thickness'].attrs['long_name'])

#using pandas rolling mean function, add a rolling mean 
rollingMeanInterval = 12 #number of months over which to take mean
tbl = pd.DataFrame({'means': PIO['mean_ice_thickness'].values})
rollingMean = np.array(tbl.rolling(rollingMeanInterval).mean()['means'])
ax.plot(PIO.time.values, rollingMean, c = 'darkblue', label = 'PIOMAS ' + str(rollingMeanInterval) + ' month mean')

#display plot with descriptive legend, ylabel, and title
plt.legend()
plt.ylabel('Mean sea ice thickness (m)')
plt.title(PIO.attrs['regions with data'] + ' mean monthly sea ice thickness')
plt.show()

## Understanding detrended uncertainty 
Detrended uncertainty captures natural variability in mean ice thickness, while ignoring the downward trend of ice thickness due to climate change. The regression line is subtracted from the mean ice thickness resulting in detrended ice thickness, which can be used to solve for natural variability (detrended uncertainty) in ice thickness.<br><br>The figure below shows these values for a single month over the period of 1978-2020 to illustrate this concept. 

Choose month of interest

month = 12 #december

Plot data

#initialize figure & axes
fig = plt.figure(figsize = [8,5]) 
ax = plt.axes()
gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) #add gridlines 
years = np.arange(1978, 1978 + len(thicknessDT[month - 1]), 1)
plt.title(PIO.attrs['regions with data'] + ': ' + pd.to_datetime(monDates[month - 1])[0].strftime('%B') + ' ice thickness & uncertainty')

#plot thickness 
monthlyThickness = PIO.sel(time = monDates[month - 1]).mean_ice_thickness.values
ax.plot(years, monthlyThickness, color = 'blue', label = 'mean ice thickness')

#plot +1 and -1 standard deviations for thickness
ax.fill_between(years, (thicknessDT[month - 1] + standDev[month - 1]) + slopes[month - 1], (thicknessDT[month - 1] - standDev[month - 1]) + slopes[month - 1], facecolor = 'dodgerblue', alpha = 0.25, edgecolor = 'none', label = 'mean uncertainty')

#plot regression line
ax.plot(years, slopes[month - 1], color = 'magenta', linestyle = '-.', label = 'regression line')

#plot detrended thickness 
ax.plot(years, thicknessDT[month - 1], color = 'maroon', label = 'mean detrended ice thickness')

#plot +1 and -1 standard deviations for detrended thickness
ax.fill_between(years, thicknessDT[month - 1] + standDevDT[month - 1], thicknessDT[month - 1] - standDevDT[month - 1], facecolor = 'salmon', alpha = 0.25, edgecolor = 'none', label = 'mean detrended uncertainty')

#add legends & show plot
plt.legend(bbox_to_anchor=(1.5, 1))
plt.show()

