# -*- coding: utf-8 -*-
""" wrangling_utils.py 
Helper functions for wrangling and analyzing data

"""

import xarray as xr 
import pandas as pd 
import numpy as np 
import numpy.ma as ma 
from scipy.interpolate import griddata


def getWinterDateRange(start_year, end_year, start_month = "November", end_month = "April"): 
    """ Gets date range for winter season/s
    Calling the function for start_year=2018, end_year=2020, start_month="November", end_month="April" will generate a date range from Nov 2018-Apr 2019 and Nov 2019-Apr 2020
    
    Args: 
        start_year (str): start year 
        end_year (str): end year 
        start_month (str, optional): month at which winter starts (default to November)
        end_month (str, optional): month at which winter ends (default to April)
        
    Returns: 
        winters (list): list of dates for all winter seasons in the input range (i.e: ['1980-11','1980-12','1981-01',
         '1981-02','1981-03','1981-04')
    """
    start_year = int(start_year)
    end_year = int(end_year)
    
    winters = []
    for year in range(start_year, end_year, 1):
        winters += pd.date_range(start = str(year) + '-' + start_month,
                                 end = str(year + 1) + '-' + end_month,
                                 freq = 'MS')
    winters = pd.to_datetime(winters)
    return winters


def is2_interp2d(is2_ds, cdr_da, method="nearest", interp_var="all"): 
    """ Perform 2D interpolation over geographic coordinates for all ICESat-2 sea ice variables with geographic coordinates in xr.Dataset
    As of 06/02/2021, xarray does not have a 2D interpolation function so this function is built on scipy.interpolate.griddata (https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html)
    This function assumes that the dataset has physical coordinates (i.e. lat,lon) as coordinates but logical coordinates (i.e. (x,y)) as dimensions (http://xarray.pydata.org/en/stable/examples/multidimensional-coords.html)
    
    Args: 
        is2_ds (xr.Dataset): ICESat-2 dataset containing variables to interpolate
        cdr_da (xr.DataArray): NSIDC sea ice concentration. Must contain the same time variable as is2_ds
        method (str,optional): interpolation method (default to "linear", choose from {‘linear’, ‘nearest’, ‘cubic’})
        interp_va (srt or list, optional): variables to interpolate (default to "all", variables with geographic coordinates)
    
    Returns: 
        ds_interp (xr.Dataset): dataset with interpolated variables
    
    """
 
    # Get 2d variables in dataset, ensure that user input was valid and raise errors if not
    def get_2d_vars(is2_ds, interp_var): 
        """ Get 2d variables in dataset for interpolation 
        """
        vars_2d = [var for var in is2_ds.data_vars if set(['latitude','longitude']).issubset(list(is2_ds[var].coords))]
        if len(vars_2d) == 0: 
            raise ValueError("The input dataset does not contain geographic coordinates that can be read. See function documentation for more information.")
        if interp_var == "all": 
            interp_var = vars_2d.copy()
        elif (type(interp_var) == str) and (interp_var in vars_2d): 
            interp_var = list(interp_var)
        else:
            raise ValueError("Invalid input for interp_var")
        return interp_var                          
    interp_var_2d = get_2d_vars(is2_ds, interp_var)
    
    # Get geographic coordinates
    lats = is2_ds['latitude'].values
    lons = is2_ds['longitude'].values
    
    # Loop through variables and timesteps and interpolate 
    np_cdr = cdr_da.values
    ds_interp = is2_ds.copy()
    for var in interp_var_2d: 
        var_interp_list = []
        da_var = is2_ds[var]
        for timestep in is2_ds.time.values: 
            da = da_var.sel(time=timestep) 
            np_cdr = cdr_da.sel(time=timestep).values
            np_da = da.values
            np_da = ma.masked_where((np.isnan(np_da)) & (np_cdr > 0.15) & (np_cdr < 1.01), np_da)
            np_interp = griddata((lons[~np_da.mask], lats[~np_da.mask]), # Interpolate
                                  np_da[~np_da.mask].flatten(),
                                  (lons, lats), 
                                  fill_value=np.nan,
                                  method=method)
            da_interp = xr.DataArray(data=np_interp, # convert numpy array --> xr.DataArray
                                     dims=da.dims, 
                                     coords=da.coords,
                                     attrs={**da.attrs,'interpolation':'interpolated from original variable using ' + method + ' interpolation'},
                                     name=da.name)
            da_interp = da_interp.where(lats < 88, np.nan) # Set pole hole to nan
            da_interp = da_interp.expand_dims("time") # Add time as a dimension. Allows for merging DataArrays 
            var_interp_list.append(da_interp)
        
        var_interp = xr.merge(var_interp_list)
        ds_interp[var] = var_interp[var] # Replace variable with interpolated variable
    
    return ds_interp


#function from regional_analysis notebook
def restrictRegionally(dataset, regionKeyList): 
    """Restrict dataset to input regions.
    
    Args: 
        dataset (xr Dataset): dataset generated by Load_IS2 notebook
        regionKeyList (list): list of region keys to restrict data to 
        
    Returns: 
        regionalDataset (xr Dataset): dataset with restricted data to input regions
    """
    
    def checkKeys(regionKeyList, regionTbl): 
        """Check that regionKeyList was defined correctly

        Raises: 
            ValueError if regionKeyList was not defined correctly 
            warning if all data was removed from the dataset
        """
        if type(regionKeyList) != list: #raise a ValueError if regionKeyList is not a list 
            raise ValueError('regionKeyList needs to be a list. \nFor example, if you want to restrict data to the Beaufort Sea, define regionKeyList = [13]')
        for key in regionKeyList: 
            if key not in list(regionTbl['key']): 
                raise ValueError('Region key ' + str(key) + ' does not exist in region mask. \n Redefine regionKeyList with key numbers from table')
        if len(regionKeyList) == 0: 
            warnings.warn('You removed all the data from the dataset. Are you sure you wanted to do this? \n If not, make sure the list regionKeyList is not empty and try again. \n If you intended to keep data from all regions, set regionKeyList = list(tbl[\"key\"])')
 
    #create a table of keys and labels
    regionMask = dataset.region_mask.attrs
    regionTbl = pd.DataFrame({'key': regionMask['keys'], 'label': regionMask['labels']})
    
    #call function to check if regionKeyList was defined correctly
    checkKeys(regionKeyList, regionTbl)
    
    #filter elements from the ice thickness DataArray where the region is the desired region
    keysToRemove = [key for key in list(regionTbl['key']) if key not in regionKeyList]
    regionalDataset = dataset.copy()
    for var in dataset.data_vars: 
        regionalVar = regionalDataset[var]
        for key in keysToRemove: 
            try:
                regionalVar = regionalVar.where(regionalVar['region_mask'] != key)
            except: 
                pass
        regionalDataset[var] = regionalVar
    
    #add new attributes describing changes made to the dataset
    labels = [regionTbl[regionTbl['key'] == key]['label'].item() for key in regionKeyList]
    if len(labels) < len(regionTbl['key']): 
        if set(regionKeyList) == set([10,11,12,13,15]): #convert to sets so unordered lists are compared
            regionalDataset.attrs['regions with data'] = 'Inner Arctic'
        else:    
            regionalDataset.attrs['regions with data'] = ('%s' % ', '.join(map(str, labels)))
        print('Regions selected: ' + regionalDataset.attrs['regions with data'])
    else: 
        regionalDataset.attrs['regions with data'] = 'All'
        print('Regions selected: All \nNo regions will be removed')
    
    return regionalDataset
