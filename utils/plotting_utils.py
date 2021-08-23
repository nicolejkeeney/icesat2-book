# +
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


# -

def plotArcticMaps(da, minval=None, maxval=None, title="", cmap="viridis", col_wrap=3): 
    """ Show data on a basemap of the Arctic. Can be one month or multiple months of data. 
    Creates an xarray facet grid. For more info, see: http://xarray.pydata.org/en/stable/user-guide/plotting.html
    
    Args: 
        da (xr DataArray): data to plot
        minval, maxval (int, optional): minimum and maximum values for the data variable (default to xarray default)
        cmap (str, optional): colormap to use (default to viridis)
        col_wrap (int, optional): number of columns of plots to display (default to 3, or None if time dimension has only one value)
    
    Returns:
        Figure displayed in notebook 
    
    """ 
    # Assign time coordinate if it doesn't exist
    try: 
        da.time
    except AttributeError: 
        da = da.assign_coords({"time":"unknown"})
    
    col_wrap = col_wrap if sum(da.time.shape) >1 else None
    col = "time" if sum(da.time.shape) >1 else None
    im = da.plot(x="longitude", y="latitude", col_wrap=col_wrap, col=col, transform=ccrs.PlateCarree(), vmin=minval, vmax=maxval, cmap=cmap, zorder=8, 
                 cbar_kwargs={'pad':0.02,'shrink': 0.8,'extend':'both'},
                 subplot_kws={'projection':ccrs.NorthPolarStereo(central_longitude=-45)})
    
    # Iterate through axes and add features 
    ax_iter = im.axes
    if type(ax_iter) != np.array: # If the data is just a single month, ax.iter returns an axis object. We need to iterate through a list or array
        ax_iter = np.array(ax_iter)
    for ax in ax_iter.flatten():
        ax.coastlines(linewidth=0.15, color = 'black', zorder = 10) # Coastlines
        ax.add_feature(cfeature.LAND, color ='0.95', zorder = 5)    # Land
        ax.add_feature(cfeature.LAKES, color = 'grey', zorder = 5)  # Lakes
        ax.gridlines(draw_labels=False, linewidth=0.25, color='gray', alpha=0.7, linestyle='--', zorder=6) # Gridlines
        ax.set_extent([-179, 179, 50, 90], crs=ccrs.PlateCarree()) # Set extent to zoom in on Arctic
       
    plt.suptitle(title, fontsize=15, horizontalalignment="center", x=0.45, y=1.04, fontweight='medium')
    plt.show()


def plotIceDrifts(uT, vT, magnitude, base_var=None, minval=None, maxval=None, col_wrap=3, res=5, scale_vec=100, vector_val=10, cmap='YlGnBu', title=""): 
    """ Plot ice drift data for one or more months on a basemap of the Arctic. 
    By default, the drift vectors will be overlayed on drift vector magnitude, but you can input a different xr.DataArray for the argument base_var to underlay a different variable. 
    
    Args: 
        uT (xr.DataArray): sea ice x velocity
        vT (xr.DataArray): sea ice y velocity
        magnitude (xr.DataArray): drift vector magnitude
        minval, maxval (int, optional): minimum and maximum values for the base variable, either base_var or magnitude if base_var is None (default to xarray default)
        base_var (xr.DataArray, optional): data to underlay below vectors (default to magnitude)
        col_wrap (int, optional): number of columns of plots to display (default to 3, or None if time dimension has only one value)
        res (int, optional): resolution of vectors (default to 5)
        scale_vec (int, optional): scaling value for displaying vectors (default to 100)
        vector_val (int, optional): value of vector to display in quiver key (default to 10 cm)
        cmap (str, optional): colormap to use (default to YlGnBu)
        title (str, optional): title to give plot (default to empty string)
    
    Returns:
        Figure displayed in notebook 
    
    """ 
    # Plot base data (will go below vectors)
    da = magnitude.copy() if base_var is None else base_var # Set base variable
    col_wrap = col_wrap if sum(da.time.shape) >1 else None
    col = "time" if sum(da.time.shape) >1 else None
    im = da.plot(x="longitude", y="latitude", col_wrap=col_wrap, col=col, transform=ccrs.PlateCarree(), vmin=minval, vmax=maxval, cmap=cmap, zorder=2, 
                 cbar_kwargs={'pad':0.02,'shrink': 0.8,'extend':'both'},
                 subplot_kws={'projection':ccrs.NorthPolarStereo(central_longitude=-45)})
    
    # Iterate through axes and add features 
    ax_iter = im.axes
    time_iter = da.time.values
    if type(ax_iter) != np.array: # If the data is just a single month, ax.iter returns an axis object. We need to iterate through a list or array
        ax_iter = np.array(ax_iter)
        timer_iter = np.array(da.time)
    
    lats = da["latitude"].values
    lons = da["longitude"].values
    for (ax,time) in zip(ax_iter.flatten(),time_iter.flatten()):
        
        # Vector computation and plotting
        u_src_crs = uT.sel(time=time).values / np.cos(lats/180*np.pi)
        v_src_crs = vT.sel(time=time).values
        magn_src_crs = np.sqrt(u_src_crs**2 + v_src_crs**2)
        var_u_scaled = u_src_crs * magnitude.sel(time=time).values / magn_src_crs
        var_v_scaled = v_src_crs * magnitude.sel(time=time).values / magn_src_crs
        Q = ax.quiver(lons[::res,::res], lats[::res,::res], var_u_scaled[::res,::res], var_v_scaled[::res,::res], 
                      transform=ccrs.PlateCarree(), units='inches', scale=scale_vec, zorder=7)
        qk = ax.quiverkey(Q, 0.85, 0.88, vector_val, str(vector_val) + ' ' + uT.attrs['units'], coordinates='axes', zorder=11, fontproperties={'size': 9})   

        # Add features and set extent 
        ax.coastlines(linewidth=0.15, color = 'black', zorder = 10) # Coastlines
        ax.add_feature(cfeature.LAND, color ='0.95', zorder = 5)    # Land
        ax.add_feature(cfeature.LAKES, color = 'grey', zorder = 5)  # Lakes
        ax.gridlines(draw_labels=False, linewidth=0.25, color='gray', alpha=0.7, linestyle='--', zorder=6) # Gridlines
        ax.set_extent([-179, 179, 50, 90], crs=ccrs.PlateCarree()) # Set extent to zoom in on Arctic

    plt.suptitle(title, fontsize=15, horizontalalignment="center", x=0.45, y=1.04, fontweight='medium')
    plt.show()


def plot_winter_means(da, start_year, end_year, da_uncertainty=None, show_uncertainty=True, title=None, colors=["blue","green","orange","red","purple","black"]):
    """ Generate a line plot of monthly means for all winters included in the dataset. DataArray must contain x,y as dimensions. 
    Plots 6 lines max. You'll need to change this function to increase the number of winters you can plot past 6 
    
    Args: 
        da (xr.DataArray): data to compute mean for 
        start_year (str): year to start computation for (i.e. "2018")
        end_year (str): year to end computation for (i.e. "2019")
        show_ucertainty (bool, optional): show detrended uncertainty (default to True)
        da_uncertainty (xr.DataArray, optional): uncertainty corresponding to da. If not inputted, the function will compute detrended uncertainty (default to None)
        title (str, optional): title to give plot (default to None)
        colors (list of str, optional): colors for each line (default to ["blue","green","orange","red","purple","black"]). Allows for max 6 lines 
    
    Returns: 
        Figure displayed in notebook 
    
    """
    
    # If user input a str for the colors input, i.e. "blue", convert to list for compatibility with function
    if type(colors) == str: 
        colors = [colors]*6
    
    # Set up plot 
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot(["Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"], np.empty((8,1))*np.nan, color=None, label=None) # Set x axis using winter months 
    gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) # Add gridlines 

    curr_yr = int(start_year)
    color_index = 0 # For looping through colors
    while curr_yr < int(end_year):
        try: 
            winter_n = pd.date_range(start="September "+str(curr_yr), end="April "+str(curr_yr+1), freq="MS") # Get winter season 
            timestamps_in_da = pd.DatetimeIndex([time for time in winter_n if time in da.time])
            winter_monthly_mean = da.sel(time=timestamps_in_da).mean(dim=["x","y"]).values # Compute mean
        except: 
            break
        
        # Plot lines 
        winter_str = "winter "+str(curr_yr)+'-'+str(curr_yr+1)
        lineplot = ax.plot(timestamps_in_da.strftime("%b"), winter_monthly_mean, marker="o", color=colors[color_index], label=winter_str)
        
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
        color_index+=1

    ax.legend(fontsize="small") # Add legend
    
    if title is not None: # Add title 
        plt.title(title)
    
    plt.show() # Display plot in notebook 


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
