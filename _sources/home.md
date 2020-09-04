ICESat-2 Sea Ice Thickness Data Analysis and Visualization
=============================================

NASA's Ice, Cloud, and Land Elevation Satellite-2 - or ICESat-2 for short - is a new satellite laser altimeter providing very high resolution elevation profiling of the entire Earth, and especially the fast-changing Polar Regions. For more information on ICESat-2, check out [https://icesat-2.gsfc.nasa.gov/](https://icesat-2.gsfc.nasa.gov/). ICESat-2 provides measurements of sea ice freeboard, the extension of sea ice above the local sea surface. Using assumptions regarding the snow depth and density ontop of the ice, and the density of the sea ice itself, these freeboard measurements can be converted to estimates of sea ice thickness [(http://www.alekpetty.com/papers/petty2020)](http://www.alekpetty.com/papers/petty2020). ICESat-2 has now collected data over two winter seasons across the entire Arctic Ocean (2018-2019 and 2019-2020) which we describe and analyze within this Jupyter Book. 


### Jupyter Book description
Jupyter Books [(jupyterbook.org)](jupyterbook.org) provide a novel means of compiling Jupyter Notebooks into one convenient and well indexed location. Jupyter Notebooks are used to provide a visual demonstration of our efforts to analyze the monthly gridded ICESat-2 winter Arctic sea ice data (sea ice freeboard and thickness), along with other relevant datasets to help us understand recent winter Arctic sea ice growth. All of the code presented in the notebooks of this book are configured to run through Google Colab, although you can download the notebooks and run them on your local drive way as well (see the [GitHub](https://github.com/nicolejkeeney/icesat2-book) repo to download the environment used).<br><br>Datasets used in the notebooks were regridded to the NSIDC North Polar Sterographic grid used by the gridded ICESat-2 thickness data and compiled into one xarray dataset, which can be simply loaded into each notebook. See the data wrangling page for more information on this process. <br><br>To run a notebook interacitvely, hover over the rocketship icon at the top of the webpage, and click Colab to open the notebook in Google Colab. 


### Author 
Nicole Keeney (NASA internship, summer 2020).
- contact: nicolekeeney@berkeley.edu
- GitHub: nicolejkeeney

### ICESat-2 contact 
Alek Petty (internship supervisor)
- contact: alek.a.petty@nasa.gov
- GitHub: akpetty


### Update history  
If you find any typos or errors in the code or have any suggestions for the book, feel free to open an issue, which you can find by mousing over the GitHub icon at the top of each page. If you are familiar with GitHub, you can also fork the book's repository and suggest an edit that way. 
 - 9/4/2020: Version 1 
