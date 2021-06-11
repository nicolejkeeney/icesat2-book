"""plotting_utils.py 
Helper functions for making plots 

"""

import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from textwrap import wrap

def plotArcticMap(da, minval = 0, maxval = 1, title = None, cbarTicks = None, cmap = 'viridis', show_mean = True, figPath = None): 
    """Plots map of the arctic on North Pole Stereo projection with one month of data overlayed, along with the sea ice edge for each month.
    Cartopy has a bug and cannot produce a contour plot on a rotated grid. Here we use a workaround from stackexchange: https://stackoverflow.com/questions/55062406/cartopy-fails-to-correctly-contour-data-on-rotated-grid
   
    Args:
        da (xr DataArray): data to plot
        minval, maxval (int): minimum and maximum values for the data variable 
        cbarTicks (list or np array of length 2): ticks to use on colorbar (default to [minval + 1, maxval +1])
        cmap (str, optional): color map (default to viridis)
        title (str, optional): title of subplots (default to name of variable)
        show_mean (bool, optional): add mean value to plot (default to True)
        figPath (str, optional): path to save fig (default to None)
        
    Returns:
        Figure displayed in notebook 
    
    """

    #initialize the figure and axes 
    proj = ccrs.NorthPolarStereo(central_longitude=-45)
    transform = ccrs.PlateCarree()
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection = proj)
    
    # colorbar settings
    cbarTicks = np.arange(minval, maxval + 1, 1) if cbarTicks is None else cbarTicks
    cbar_label = "" #create colorbar label using variable name and units if available as attributes
    attrs = da.attrs
    if 'long_name' in attrs:
        cbar_label += attrs['long_name'] + " "
    if 'units' in attrs:
        cbar_label += "(" + attrs['units'] + ")"
    elif 'unit' in attrs: 
        cbar_label += "(" + attrs['unit'] + ")"

    #plot the data
    da.plot(x = 'longitude', y = 'latitude', vmin = minval, vmax = maxval, extend = 'both', 
            ax = ax, add_colorbar = True, transform = transform, zorder = 2, cmap = cmap, 
            cbar_kwargs = {'label': "\n".join(wrap(cbar_label, 50)), 'orientation': 'horizontal', 'shrink': 0.75, 'pad': 0.025})
    
    #add features to the map
    ax.coastlines(linewidth=0.15, color = 'black', zorder = 10) #add coastlines 
    ax.add_feature(cfeature.LAND, color ='0.95', zorder = 5) #add land 
    ax.add_feature(cfeature.LAKES, color = 'grey', zorder = 5) #add lakes 
    ax.gridlines(draw_labels = False, linewidth = 0.25, color = 'gray', alpha = 0.7, linestyle = '--', zorder = 6) #add gridlines
    ax.set_extent([-179, 179, 55, 90], crs = transform) #zoom in so map only displays the Arctic
    
    if title != None: #add title
        ax.set_title("\n".join(wrap(title, 40)), fontsize = 13, y = 1, fontweight = 'medium')
    
    if show_mean: #add mean value 
        mean_val = round(da.mean().values.item(),2)
        ax.text(x = 0.68, y = 0.97, s = "mean val: " + str(mean_val), zorder = 20, horizontalalignment = 'left', verticalalignment = 'top', fontsize = 10, transform = ax.transAxes, bbox = {'facecolor':'white', 'edgecolor':'lightgrey', 'boxstyle':'round'})
    
    if figPath != None: # save and display figure
        plt.savefig(figPath, dpi = 200)
    plt.show()


def polar_fig_and_axes(hemisphere, figsize=(10,10), num_rows=1, num_columns=3, lat=55): 
    """ Generate empty figure and axes for overlaying data 
    
    Args: 
        hemisphere ('nh' or 'sh'): string indicating projection to use for map
        figsize (tuple, optional): size to use for figure (default to (10,10))
        num_rows (int, optional): number of axes to generate (default to 3)
        num_columns (int, optional): number of columns to generate (default to 1)
        lat (float, optional): positive float value to restrict above (northern hemisphere) or below (southern hemisphere) (default to 55)
    
    Returns: 
        fig (matplotlib.figure.Figure): figure
        axes (tuple of cartopy.mpl.geoaxes.GeoAxesSubplot): axis for each subplot
    """
    # Check if inputs are valid
    if any(type(val) not in [float, int] for val in [num_rows, num_columns, lat]): 
        raise ValueError("num_rows, num_columns, and lat need to be a valid float or integer value.")
    
    if lat < 0: # Lat needs to be a positive number for the southern hemisphere; negative value already coded in (see code chunk below for defining extent)
        lat = lat*(-1)
    
    if lat > 90: 
        raise ValueError("Input a valid value for latitude in [0,90]")
    
    # Define projection and extent depending on hemisphere of interest
    if hemisphere == 'nh': 
        projection = ccrs.NorthPolarStereo()
        extent = [-179,179,lat,90]
    elif hemisphere == 'sh': 
        projection = ccrs.SouthPolarStereo()
        extent = [-180,180,-90,-(lat)]
        
    # Initialize figure and subplot axes 
    fig, axes = plt.subplots(num_rows, num_columns, figsize = figsize, subplot_kw = {'projection': projection})
    
    try: 
    # Flatten axes because if num_column > 1, axes is a 2D array that cannot be iterated through
        axes = axes.flatten() 
    except: 
        # If only one axis is generated, convert axis object to an array so it can be iterated through 
        axes = np.array([axes])

    for ax in axes: # Make axes pretty 
        ax.coastlines(linewidth = 0.25, color = 'black', zorder = 10) #add coastlines 
        ax.set_extent(extent, ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color ='0.95', zorder = 5) #add land 
        ax.add_feature(cfeature.LAKES, color = 'grey', zorder = 5) #add lakes 
        ax.gridlines(draw_labels = False, linewidth = 0.25, color = 'grey', alpha = 0.75, linestyle='--', zorder = 15)
    
    return fig, axes


def arcticComparisonMaps(data1, data2, plot_diff=True, vmin=None, vmax=None, vmin_diff=None, vmax_diff=None, hemisphere='nh', cmap="viridis", lat=55, title1="data_1", title2="data_2"):
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
        hemisphere ('nh' or 'sh', optional): string indicating projection to use for map (default to 'nh')
        title1 (str, optional): title to give plot of data1 (default to "data1")
        title2 (str, optional): title to give plot of data2 (default to "data2")
        cmap (str, optional): colormap to use (default to viridis)
        lat (float, optional): positive float value to restrict above (northern hemisphere) or below (southern hemisphere) (default to 55)
    
    Returns: 
        figure displayed in notebook 
        
    """

    # Compute min and max for plotting 
    if vmin is None: 
        vmin = round(np.nanpercentile(data1.values, 1),1)
    if vmax is None: 
        vmax = round(np.nanpercentile(data1.values, 99),1) 
        
    # Generate figure and axes 
    num_columns = 3 if (plot_diff==True) else 2
    fig, axes = polar_fig_and_axes(hemisphere, figsize=(16,6), num_columns=num_columns, lat=lat)
    titles = [title1,title2]

    # Overlay data
    im1 = data1.plot(x='longitude', y='latitude', ax=axes[0], cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=8, cbar_kwargs = {'orientation':'horizontal', 'pad':0.02, 'extend':'both'})
    im2 = data2.plot(x='longitude', y='latitude', ax=axes[1], cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=8, cbar_kwargs = {'orientation':'horizontal', 'pad':0.02, 'extend':'both'})
    
    # Plot gricell difference 
    if plot_diff: 
        gridcell_diff = data2-data1
        if (vmax_diff is not None) and (vmin_diff is not None): 
            diff_im = gridcell_diff.plot(x='longitude', y='latitude', ax=axes[2], transform=ccrs.PlateCarree(), vmin = vmin_diff, vmax=vmax_diff, zorder=8, cmap='coolwarm', cbar_kwargs={'orientation':'horizontal', 'pad':0.02, 'extend':'both', 'label':'difference'})
        else: 
            diff_im = gridcell_diff.plot(x='longitude', y='latitude', ax=axes[2], transform=ccrs.PlateCarree(), zorder=8, center=0, cmap='coolwarm', cbar_kwargs={'orientation':'horizontal', 'pad':0.02, 'extend':'both', 'label':'difference'})
        titles.append("Griddcell Difference")
    
    for ax, title in zip(axes,titles): # Add titles to figure and axes 
        ax.set_title("\n".join(wrap(title, 40)), fontsize = 13, y = 1, fontweight = 'medium')

    plt.show() # Show fig in notebook 
