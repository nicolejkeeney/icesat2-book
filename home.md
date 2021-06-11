ICESat-2 Sea Ice Thickness Data Analysis and Visualization
=============================================

NASA's Ice, Cloud, and Land Elevation Satellite-2 (ICESat-2) is a new satellite laser altimetry mission providing very high resolution elevation profiling of the entire Earth, and especially the fast-changing Polar Regions. For more information on ICESat-2, check out [https://icesat-2.gsfc.nasa.gov/](https://icesat-2.gsfc.nasa.gov/). ICESat-2 provides measurements of sea ice freeboard, the extension of sea ice above the local sea surface. Using assumptions regarding the snow depth and density ontop of the ice, and the density of the sea ice itself, these freeboard measurements can be converted to estimates of sea ice thickness [(http://www.alekpetty.com/papers/petty2020)](http://www.alekpetty.com/papers/petty2020). ICESat-2 has now collected data over two winter seasons across the entire Arctic Ocean (2018-2019 and 2019-2020) which we describe and analyze within this Jupyter Book. 

<!-- #region -->
## Jupyter Book description
Jupyter Books [(jupyterbook.org)](https://jupyterbook.org/intro.html) provide a novel means of compiling Jupyter Notebooks into one convenient and well-indexed location. Jupyter Notebooks are used to provide a visual demonstration of our efforts to analyze the monthly gridded ICESat-2 winter Arctic sea ice data (sea ice freeboard and thickness), along with other relevant datasets to help us understand recent winter Arctic sea ice growth. All of the code presented in the notebooks of this book are configured to run in Binder, although you can download the notebooks and run them on your local drive way as well (see the [GitHub](https://github.com/nicolejkeeney/icesat2-book) repo to download the environment used). Using Binder you can run the code without needing to download anything-- even python! To run a notebook in Binder, just click Binder under the rocketship icon at the top of each notebook.<br><br>

The ICESat-2 monthly gridded data can be downloaded from the google storage bucket associated with this jupyter book, along with additional datasets used in the notebook. All other datasets used were regridded to the NSIDC North Polar Sterographic grid (the grid used by the ICESat-2 data) such that they can be easily compared between each other. See the data wrangling page for more information on this process.<br><br> 

This book was created by the author during a summer 2020 student internship at NASA Goddard, and updated in June 2021 for the publication of Petty et al. (2021). 

## A note on xarray 
All of the notebooks in this notebook utilize [xarray](http://xarray.pydata.org/en/stable/), a python package built for working with multi-dimensional data like the monthly gridded sea ice data. Xarray is especially useful for time series data and allows for easily plotting data on map projections via compatability with the python package cartopy. 


## Author 
Nicole Keeney
- contact: nicolekeeney@berkeley.edu
- GitHub: nicolejkeeney

## ICESat-2 contact 
Alek Petty
- contact: alek.a.petty@nasa.gov
- GitHub: akpetty


## Update history  
If you find any typos or errors in the code or have any suggestions for the book, feel free to open an issue, which you can find by mousing over the GitHub icon at the top of each page. If you are familiar with GitHub, you can also fork the book's repository and suggest an edit that way. 
 - 9/4/2020: Version 1 
 - 11/18/2020: Updated with new gridded ice thickness dataset
<!-- #endregion -->

```python

```
