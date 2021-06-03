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


def compute_stats(data): 
    """ Compute mean value and standard deviation for each time step
    If data is a xr.Dataset containing more than one data variable, the funciton will convert it to a xr.DataArray and compute statistics using the first data variable
    
    Args: 
        data (xr.Dataset or xr.DataArray)
    
    Returns: 
        mean (xr.DataArray): mean value 
        std (xr.DataArray): standard deviation 
    """
    
    # Turn Dataset object into DataArray
    try: 
        data = data[list(data.data_vars)[0]]  # Get first data variable from list of data variables 
    except: # If data is already a DataArray, pass
        pass
    
    # Compute mean value
    mean = data.mean(dim = ['lat','lon'], skipna = True)
    mean = mean.assign_attrs({'long_name':'mean value', 'reference': 'http://xarray.pydata.org/en/stable/generated/xarray.DataArray.mean.html'})
    mean.name = 'mean_value'

    # Compute standard deviation
    std = data.std(dim = ['lat','lon'], skipna = True)
    std = std.assign_attrs({'long_name':'standard deviation', 'reference': 'http://xarray.pydata.org/en/stable/generated/xarray.DataArray.std.html'})
    std.name = 'std'

    return mean, std



def compute_comparative_stats(obs_data, model_data): 
    """ Compute comparative statistics 
    Requires that both input datasets are on the same grid and contain the dimensions "time", "lat", "lon"
    Statistical metrics computed using xskillscore package.
    If obs_data or model_data is a xr.Dataset containing more than one data variable, the funciton will convert it to a xr.DataArray and compute statistics using the first data variable
    
    Args: 
        obs_data (xr.Dataset or xr.DataArray): observational data 
        model_data (xr.Dataset or xr.DataArray): model data 
        
    Returns: 
        gridcell_diff (xr.DataArray): gridcell difference model - observational data 
        RMSE_monthly (xr.DataArray): root mean square error per month; returns one value per month 
        RMSE_tot (xr.DataArray): root mean square error over all dimensions; returns a single value 
        pearson_monthly (xr.DataArray): pearson r per month; returns one value per month
        pearson_tot (xr.DataArray): pearson r over all dimensions; returns a single value
    
    """
    # Turn Dataset objects into DataArrays
    try: 
        model_data = model_data[list(model_data.data_vars)[0]]  # Get first data variable from list of data variables 
    except: # If data is already a DataArray, pass
        pass
    try: 
        obs_data = obs_data[list(obs_data.data_vars)[0]]  # Get first data variable from list of data variables 
    except: # If data is already a DataArray, pass
        pass
    
    # Compute griddcell_diff
    gridcell_diff = (model_data - obs_data).assign_attrs({'long_name': 'gridcell difference (model - observation)'})
    gridcell_diff.name = 'gridcell_diff'
    
    # Rechunk to avoid chunking issues 
    try: 
        obs_data = obs_data.chunk({"time":-1})
        model_data = model_data.chunk({"time":-1})
    except: 
        pass 
    
    # Compute RMSE 
    rmse_tot = xs.rmse(obs_data, model_data, dim = ["time","lat", "lon"], skipna = True).assign_attrs({'long_name':'root mean square error', 'reference':'https://xskillscore.readthedocs.io/en/stable/api/xskillscore.rmse.html#xskillscore.rmse'})
    rmse_monthly = xs.rmse(obs_data, model_data, dim = ["lat", "lon"], skipna = True).assign_attrs({'long_name':'root mean square error', 'reference':'https://xskillscore.readthedocs.io/en/stable/api/xskillscore.rmse.html#xskillscore.rmse'})
    rmse_tot.name = 'RMSE'
    rmse_monthly.name = 'RMSE_monthly'
    
    # Compute pearson R
    pearson_tot = xs.pearson_r(obs_data, model_data, dim = ["time", "lat", "lon"], skipna = True).assign_attrs({'long_name':'pearson r', 'reference':'https://xskillscore.readthedocs.io/en/stable/api/xskillscore.pearson_r.html#xskillscore.pearson_r'})
    pearson_monthly = xs.pearson_r(obs_data, model_data, dim = ["lat", "lon"], skipna = True).assign_attrs({'long_name':'pearson r', 'reference':'https://xskillscore.readthedocs.io/en/stable/api/xskillscore.pearson_r.html#xskillscore.pearson_r'})
    pearson_tot.name = 'pearson_R'
    pearson_monthly.name = 'pearson_R_monthly'
    
    return gridcell_diff, rmse_tot, rmse_monthly, pearson_tot, pearson_monthly