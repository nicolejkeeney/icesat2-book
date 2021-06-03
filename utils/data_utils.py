""" data_utils.py 
Helper functions for reading and downloading data 

"""

import xarray as xr 
import os 
import pandas as pd


def read_is2_data(data_dir="IS2SITMOGR4"): 
    """ Read in ATLAS/ICESat-2 Monthly Gridded Sea Ice Freeboard dataset. 
    If the file does not already exist on the user's local drive, it is downloaded from the books google storage bucket (https://console.cloud.google.com/storage/browser/is2-pso-seaice)
    The netcdf files for each month are then read in as an xr.Dataset object
    
    Args: 
        data_dir (str, optional): name of data directory containing ICESat-2 data (default to "IS2SITMOGR4", the name of the directory in the bucket)
    Returns: 
        is2_ds (xr.Dataset): data 
    
    """
    # Download data from bucket if it doesn't exist on the user's local drive
    exists_locally = os.path.isdir(data_dir)
    if (exists_locally == False): 
        print("Downloading ICESat-2 data from the google storage bucket...")
        os.system("gsutil -m cp -r gs://is2-pso-seaice/" + data_dir + " ./") # Make sure theres a space before the final ./ (i.e. " ./")
        print("Download complete")

    # Read in files for each month as a single xr.Dataset
    # Need to create a preprocessing function to call before merging because dimensions and coordinates are not set
    # This allows each DataArray for each month to be merged into one xr.Dataset
    def xr_set_coords_and_dims(da_monthly):
        da_monthly = da_monthly.set_coords(["latitude","longitude","xgrid","ygrid"]) # Set data variables as coordinates
        da_monthly = da_monthly.expand_dims("time") # Set month as a dimension 
        return da_monthly
    
    filenames = os.listdir(data_dir)
    is2_ds = xr.open_mfdataset([data_dir + "/" + filename for filename in filenames], # Filepath, including data directory in path
                               concat_dim=["time"], 
                               combine='nested', 
                               preprocess=xr_set_coords_and_dims)
    time = [file.split("IS2SITMOGR4_01_")[1].split("_004_001.nc")[0] for file in filenames] # Get time from filenames
    is2_ds = is2_ds.assign_coords({"time":pd.to_datetime(time, format = "%Y%m")}) # Add time as coordinate
    return is2_ds


def read_book_data(filepath = "icesat2-book-data.nc"): 
    """ Read in data for ICESat2 jupyter book. 
    If the file does not already exist on the user's local drive, it is downloaded from the books google storage bucket (https://console.cloud.google.com/storage/browser/is2-pso-seaice)
    The netcdf file is then read in as an xr.Dataset object 
    
    Args: 
        filepath (str, optional): name of file to read in (default to "icesat2-book-data.nc", the name of the file in the bucket)
    Returns: 
        book_ds (xr.Dataset): data 
    
    """
    exists_locally = os.path.isfile(filepath) # Check if file exists on local drive
    if (exists_locally == False): # Download data 
        print("Downloading jupyter book data from the google storage bucket...")
        os.system("gsutil -m cp gs://is2-pso-seaice/" + filepath + " ./") # Make sure theres a space before the final ./ (i.e. " ./")
        print("Download complete")

    book_ds = xr.open_dataset(filepath)
    return book_ds


