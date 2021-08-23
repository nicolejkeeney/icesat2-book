# Plotting

This is a markdown rendering of the plotting_utils module used in the notebooks. It is provided here for user reference. The `plotting_utils.py` module can be viewed and downloaded from the github repository. 


```
""" plotting_utils.py 

    Functions that create nice plots to visualize the data 

"""

import numpy as np 
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from textwrap import wrap
import matplotlib as mpl
import matplotlib.pyplot as plt
```


```
def plot_winter_means(da, start_year, end_year, da_uncertainty=None, show_uncertainty=True, title=None):
    """ Generate a line plot of monthly means for all winters included in the dataset. DataArray must contain x,y as dimensions. 
    
    Args: 
        da (xr.DataArray): data to compute mean for 
        start_year (str): year to start computation for (i.e. "2018")
        end_year (str): year to end computation for (i.e. "2019")
        show_ucertainty (bool, optional): show detrended uncertainty (default to True)
        da_uncertainty (xr.DataArray, optional): uncertainty corresponding to da. If not inputted, the function will compute detrended uncertainty (default to None)
        title (str, optional): title to give plot (default to None)
        
    Returns: 
        Figure displayed in notebook 
    
    """
    
    # Set up plot 
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot(["Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"], np.empty((8,1))*np.nan, color=None, label=None) # Set x axis using winter months 
    gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) # Add gridlines 

    curr_yr = int(start_year)
    while curr_yr < int(end_year):
        try: 
            winter_n = pd.date_range(start="September "+str(curr_yr), end="April "+str(curr_yr+1), freq="MS") # Get winter season 
            timestamps_in_da = pd.DatetimeIndex([time for time in winter_n if time in da.time])
            winter_monthly_mean = da.sel(time=timestamps_in_da).mean(dim=["x","y"]).values # Compute mean
        except: 
            break
        
        # Plot lines 
        winter_str = "winter "+str(curr_yr)+'-'+str(curr_yr+1)
        lineplot = ax.plot(timestamps_in_da.strftime("%b"), winter_monthly_mean, marker="o", label=winter_str)
        
        if show_uncertainty: 
            if da_uncertainty is None: 
                da_mean_unc = da.sel(time=timestamps_in_da).std(dim=["x","y"]).values
                unc_descrip = '+/-1 std'
            else: 
                da_mean_unc = da_uncertainty.sel(time=timestamps_in_da).mean(dim=["x","y"]).values
                unc_descrip = 'uncertainty'
            ax.fill_between(timestamps_in_da.strftime("%b"), winter_monthly_mean-da_mean_unc, winter_monthly_mean+da_mean_unc, alpha=0.1, edgecolor=None, label=winter_str+" "+unc_descrip, color=lineplot[0].get_color())  

        try: # Try adding attributes for name and units as ylabel
            ylabel_str = "\n".join(wrap(da.attrs["long_name"]+" ("+da.attrs["units"]+")", 40))
            plt.ylabel(ylabel_str)
        except: # Leave ylabel blank if not 
            pass
        
        curr_yr+=1

    ax.legend(fontsize="small") # Add legend
    
    if title is not None: # Add title 
        plt.title(title)
    
    plt.show() # Display plot in notebook 
```


```

```


```
def arcticComparisonMaps(data1, data2, plot_diff=True, vmin=None, vmax=None, vmin_diff=None, vmax_diff=None, cmap="viridis", title1="data_1", title2="data_2"):
    """ Plot comparison maps for two xr.DataArrays, along with gridcell difference 
    Both input DataArrays must be on the same grid and contain the coordinates latitude and longitude
    
    Args: 
        data1 (xr.DataArray): data to plot on leftmost plot. 
        data2 (xr.DataArray): data to plot on middle plot
        plot_diff (bool, optional): plot gridcell difference (data1 - data2) on a third map? (default to True)
        vmin (float, optional): minimum value to use colorbar (default to 1st percentile of data1)
        vmax (float, optional): maximum value to use colorbar (default to 99th percentile of data1)
        vmin_diff (float, optional): minimum value to use for plotting data range (default to mapping default, with 0 at center)
        vmax_diff (float, optional): maximum value to use for plotting data range (default to mapping default, with 0 at center)
        title1 (str, optional): title to give plot of data1 (default to "data1")
        title2 (str, optional): title to give plot of data2 (default to "data2")
        cmap (str, optional): colormap to use (default to viridis)
        
    Returns: 
        figure displayed in notebook 
        
    """

    # Compute min and max for plotting 
    vmin = vmin if vmin is not None else round(np.nanpercentile(data1.values, 1),1)
    vmax = vmax if vmax is not None else round(np.nanpercentile(data1.values, 99),1) 
        
    # Plot
    titles = [title1,title2]
    num_columns = 3 if (plot_diff==True) else 2
    fig, axes = plt.subplots(1, num_columns, figsize=(16,6), subplot_kw={'projection':ccrs.NorthPolarStereo(central_longitude=-45)})
    im1 = data1.plot(x='longitude', y='latitude', ax=axes[0], cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=8, cbar_kwargs = {'orientation':'horizontal', 'pad':0.02, 'extend':'both'})
    im2 = data2.plot(x='longitude', y='latitude', ax=axes[1], cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=8, cbar_kwargs = {'orientation':'horizontal', 'pad':0.02, 'extend':'both'})
     
    if plot_diff: # Plot gricell difference
        gridcell_diff = data2 - data1
        diff_im = gridcell_diff.plot(x='longitude', y='latitude', ax=axes[2], transform=ccrs.PlateCarree(), vmin=vmin_diff, vmax=vmax_diff, zorder=8, cmap='coolwarm', cbar_kwargs={'orientation':'horizontal','pad':0.02,'extend':'both','label':'difference'})
        titles.append("Griddcell Difference")
    
    for ax, title in zip(axes,titles): # Add title and features to axes 
        ax.set_title("\n".join(wrap(title, 40)), fontsize = 13, y = 1, fontweight = 'medium')
        ax.coastlines(linewidth=0.15, color = 'black', zorder = 10) # Coastlines
        ax.add_feature(cfeature.LAND, color ='0.95', zorder = 5)    # Land
        ax.add_feature(cfeature.LAKES, color = 'grey', zorder = 5)  # Lakes
        ax.gridlines(draw_labels=False, linewidth=0.25, color='gray', alpha=0.7, linestyle='--', zorder=6) # Gridlines
        ax.set_extent([-179, 179, 50, 90], crs=ccrs.PlateCarree()) # Set extent to zoom in on Arctic

    plt.show() # Show fig in notebook 
```
