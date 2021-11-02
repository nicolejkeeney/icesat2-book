ICESat-2 Sea Ice Thickness Data Analysis and Visualization
=============================================

NASA's Ice, Cloud, and Land Elevation Satellite-2 (ICESat-2) is a new satellite laser altimetry mission providing high resolution elevation profiling of the entire Earth's surface, particularly in the fast-changing Polar Regions. ICESat-2 provides measurements of sea ice freeboard, the extension of sea ice above the local sea surface. Using assumptions regarding the snow depth and density on top of the ice, along with the density of the sea ice itself, these freeboard measurements can be converted to estimates of sea ice thickness [(http://www.alekpetty.com/papers/petty2020)](http://www.alekpetty.com/papers/petty2020). Since its launch in 2018, ICESat-2 has collected and released data over several winter seasons across the entire Arctic Ocean, which we describe and analyze within this Jupyter Book. <br><br> For more information on ICESat-2, see the project homepage: [https://icesat-2.gsfc.nasa.gov/](https://icesat-2.gsfc.nasa.gov/).


# Jupyter Book description
[Jupyter Books](https://jupyterbook.org/intro.html) provide a novel means of compiling Jupyter notebooks into one convenient and well-indexed location. Here, jupyter notebooks are used to provide a visual demonstration of our efforts to analyze the monthly gridded ICESat-2 winter Arctic sea ice data (freeboard and thickness), along with other relevant datasets to help us understand recent winter Arctic sea ice growth.<br><br>We've also set up the book so that users can easily run the code without needing to download anything-- even python-- by using a hosting service called Binder. To run a notebook (chapter pages in the book) in Binder, just click the **Binder** tab under the rocket ship icon at the top of each notebook. This option is configured for all notebooks except the modules in the Helper Functions section and the data wrangling notebook. 


# Contact 

**Nicole Keeney (author)**
<br>nicolejkeeney@gmail.com
<br>GitHub: nicolejkeeney

**Alek Petty (ICESat-2 contact)**
<br>alek.a.petty@nasa.gov
<br>GitHub: akpetty

<!-- #region -->
# Accessing the data 

ICESat-2 monthly gridded data can be downloaded from the [google storage bucket](https://console.cloud.google.com/storage/browser/sea-ice-thickness-data) associated with this Jupyter Book as individual netcdf files for each month. We've also generated a single netcdf file, stored in the same bucket under the name `icesat2-book-data.nc`, that contains all the monthly data. The file also contains all the other datasets used in the notebook as data variables, which we use in the book to help contextualize the sea ice and atmospheric conditions of each winter. All datasets included used were regridded to the NSIDC North Polar Sterographic grid (the native grid of the ICESat-2 data used), such that they can be easily compared with each other. See the data wrangling page for more information each dataset and on on the regridding process.<br><br> 


ICESat-2 data is also publicly avaiable through the [National Snow and Ice Data Center (NSIDC)](https://nsidc.org/data/icesat-2)
<!-- #endregion -->

# Update history  
If you find any typos or errors in the code or have any suggestions for the book, feel free to open an issue, which you can find by mousing over the GitHub icon at the top of each page. If you are familiar with GitHub, you can also fork the book's repository and suggest an edit that way. 
 - 9/4/2020: Version 1
 - 11/18/2020: Updated with version 2 ICESat-2 data product for [AGU Fall 2020 poster highlighting the book](https://agu.confex.com/agu/fm20/meetingapp.cgi/Paper/684153). 
 - 6/14/2021: Updated for the publication of Petty et al. (2021) with final ICESat-2 data product and cleaner notebooks. Transitioned from Google Colab interactivity to Binder. 
 - 10/25/2021: Added interactive plotting using hvplot. Improved interpolation/smoothing method for ICESat-2 data and added notebook to demonstrate steps. 


# A note on xarray 
All of the notebooks in this notebook utilize [xarray](http://xarray.pydata.org/en/stable/), a python package built for working with multi-dimensional data like the monthly gridded sea ice data. Xarray is especially useful for time series data and allows for easily plotting data on map projections via compatability with the python packages cartopy and hvplot. 
