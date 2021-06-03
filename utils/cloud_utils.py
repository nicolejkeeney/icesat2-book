"""cloud_utils.py 

    Functions for reading zarr files from a google cloud bucket. 
    
    Author: Nicole Keeney 
    Date Created: 03-24-2021
    Modification History: 
        - Removed project_name input since it's not required to interact with the bucket. Added ability to upload files to a public bucket that doesn't require a key (05-25-2021)

"""

import os
import sys
import xarray as xr 
import pandas as pd
import gcsfs 


def create_file_system(token_path = None):
    """ Create GCSFileSystem. 
    Used to read and write zarr files to the google bucket.
    Uses gcsfs package for python interface with google cloud.
    
    Args: 
        token_path (str): path to JSON service key on local drive
    
    Returns: 
        fs (GCSFileSystem)
        
    Raises: 
        ValueError
        FileNotFoundError
    """
    
    try: 
        fs = gcsfs.GCSFileSystem(token = token_path)
        return fs
    except ValueError as err: 
        err.args = ("File " + token_path + " is not a path to a valid JSON file.",)
        raise
    except FileNotFoundError as err: 
        err.args = ("File '" + token_path + "' not found.",)
        raise
        
        
def create_mutable_mapping(bucket_name, data_path, fs, write = False): 
    """ Create mutable mapping object.
    Used to read and write zarr files to the google bucket.
    Uses gcsfs package for python interface with google cloud.
    Set write = True to write a file to the bucket, and write = False to read a file from the bucket 
    
    Args: 
        bucket_name (str): name of google bucket 
        data_path (str): path to the directory or dataset in the bucket 
        fs (GCSFileSystem): file system object 
        write (bool, optional): read (False) or write (True) the input data_path to the bucket? (default to False)
        
    
    Returns: 
        gcsmap (GCSMap)
    
    """
    # Read a file from the bucket
    if write == False: 
        try: 
            gcsmap = gcsfs.mapping.GCSMap(bucket_name + '/' + data_path, gcs = fs, check = True, create = False)
        except ValueError as err: 
            err.args = ("Path " + data_path + " does not exist.",)
            raise 
    
    # Write a file to the bucket
    elif write == True:
        gcsmap = gcsfs.mapping.GCSMap(bucket_name + '/' + data_path, gcs = fs, check = False, create = True)
    
    return gcsmap
        
        
        
def list_files_in_bucket(bucket_name, data_path, fs): 
    """ List all the files in the google bucket at an input data path 
    
    Args: 
        bucket_name (str): name of google bucket 
        data_path (str): path to the directory or dataset in the bucket 
        fs (GCSFileSystem): file system object 
    
    Returns: 
        files_in_path (list): list of filenames 
    """
    
    try: 
        files_in_path = fs.ls(bucket_name + "/" + data_path)
        return files_in_path
    
    except FileNotFoundError as err: 
        err.args = ("The input filepath '" + data_path + "' does not exist, or the bucket name '" + bucket_name + "'is wrong.",)
        raise 
        
        
def upload_ds_to_gc(xr_ds, zarr_filename, bucket_name, cloud_dir = "", token_path = None): 
    """ Upload an xr.Dataset to google cloud bucket as a zarr
    
    Args: 
        xr_ds (xr.Dataset): dataset to convert to zarr and upload to the bucket
        cloud_dir (str): location on the cloud where you want to upload data. Must be a directory
        zarr_filename (str): name you want the file to have on google drive
        bucket_name (str): name of bucket on google cloud
        token_path (str, optional): path to token on local drive. This allows you to interact with a private bucket (default to None)
    
    Returns: 
        Writes xr_ds to google bucket
        
    """
    
    # ------------- Check that inputs are in the correct format -------------
    
    # Add slash to directory if not already there 
    if len(cloud_dir) != 0: 
        if cloud_dir[-1] != "/": 
            cloud_dir = cloud_dir + '/'
            print("Cloud directory path changed to " + cloud_dir)
    
    # Make sure zarr filename includes zarr as the extension. 
    # If not, add zarr as the extension
    extension = os.path.splitext(zarr_filename)[1]
    if extension != ".zarr":
        zarr_filename = zarr_filename + ".zarr"
        print("Zarr filename changed to " + zarr_filename)

    if type(xr_ds) == xr.core.dataarray.DataArray: 
        xr_ds = xr_ds.to_dataset()
        
    # ------------- Set up cloud stuff -------------
    
    # Create file system
    fs = create_file_system(token_path = token_path)
    
    # Check that the bucket exists 
    files_in_bucket = list_files_in_bucket(bucket_name = bucket_name, data_path = cloud_dir, fs = fs)

    # Check that a file by the same name doesn't already exist in the google bucket 
    if any([s for s in files_in_bucket if zarr_filename in s]): 
        print("WARNING: A file by the name " + zarr_filename + " already exists in the bucket. \nThis file will replace the already existing file by the same name.")

    # Create zarr folder in cloud (create = True)
    gcsmap = create_mutable_mapping(bucket_name = bucket_name, data_path = cloud_dir + zarr_filename, fs = fs, write = True)
   

    # ------------- Write file to the bucket -------------
    
    xr_ds.to_zarr(store = gcsmap, consolidated = True, mode = 'w')         
    print("File written to " + cloud_dir + zarr_filename)
    
    
                  

def get_gc_zarr(bucket_name, data_path, token_path = None, chunks = "auto"): 
    """ Get zarr from google cloud bucket 
    
    Args: 
        bucket_name (str): name of bucket within google cloud project, containing data to be accessed 
        token_path (str): path to JSON service key on local drive
        data_path (str): path to data on bucket, not including bucket name in the path 
        chunks (int, dict, tuple, None, or str): chunk sizes along each dimension (default to "auto") 
    
    Returns: 
        data (xarray dataset): data chunked according to chunks argument
    
    Raises: 
        FileNotFoundError: Input token path does not exist in local directory
        JSONDecodeError: JSON file is not valid
        ValueError: Data does not exist in bucket in input path (data_path)
        
    """
    # Create FileSystem
    fs = create_file_system(token_path = token_path)
    
    # Create mutable mapping 
    gcsmap = create_mutable_mapping(bucket_name = bucket_name, data_path = data_path, fs = fs, write = False)
    
    # Load data from bucket
    data = xr.open_zarr(store = gcsmap, chunks = chunks) # Open zarr
    return data