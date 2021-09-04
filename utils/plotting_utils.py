""" plotting_utils.py 

    Functions that create nice plots to visualize the data 

"""
import hvplot.xarray
import hvplot.pandas
import numpy as np 
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from textwrap import wrap
import matplotlib as mpl
import matplotlib.pyplot as plt


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


def polar_fig_and_axes(hemisphere, figsize = (10,10), num_rows = 1, num_columns = 3, lat = 55, show_map_features=True): 
    """ Generate empty figure and axes for overlaying data 
    
    Args: 
        hemisphere ('nh' or 'sh'): string indicating projection to use for map
        figsize (tuple, optional): size to use for figure (default to (10,10))
        num_rows (int, optional): number of axes to generate (default to 3)
        num_columns (int, optional): number of columns to generate (default to 1)
        lat (float, optional): positive float value to restrict above (northern hemisphere) or below (southern hemisphere) (default to 55)
        show_map_features (bool, optional): show land on map? (default to True)
    
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
        projection = ccrs.NorthPolarStereo(central_longitude=-45)
        extent = [-179,179, lat, 90]
    elif hemisphere == 'sh': 
        projection = ccrs.SouthPolarStereo()
        extent = [-180,180, -90, -(lat)]
        
    # Initialize figure and subplot axes 
    fig, axes = plt.subplots(num_rows, num_columns, figsize=figsize, subplot_kw={'projection':projection})
    
    try: 
    # Flatten axes because if num_column > 1, axes is a 2D array that cannot be iterated through
        axes = axes.flatten() 
    except: 
        # If only one axis is generated, convert axis object to an array so it can be iterated through 
        axes = np.array([axes])

    # Make axes pretty 
    for ax in axes: 
        ax.set_extent(extent, ccrs.PlateCarree())
        ax.gridlines(draw_labels = False, linewidth=0.25, color='grey', alpha=0.75, linestyle='--', zorder=15)
        if show_map_features == True: 
            ax.coastlines(linewidth = 0.25, color='black', zorder=10) #add coastlines 
            ax.add_feature(cfeature.LAND, color='0.95', zorder=5) #add land 
            ax.add_feature(cfeature.LAKES, color='grey', zorder=5) #add lakes 
    
    return fig, axes


def compute_vmin_vmax(da): 
    vmin = np.nanpercentile(da.values, 1)
    vmax = np.nanpercentile(da.values, 99)
    return vmin, vmax


# +
# These functions are all used in the polarComparisonMaps function. Its a beefy function, meant to handle a lot of use cases! 

def add_time_dim(da): 
    """ Add time as a dimension 
    xarray is annoying with how it treats sel vs. isel when indexing items.
    isel returns time as a coord, but not a dimension. sel returns time as both a coord and dimension. 
    This code chunk account for that, and sets time as a dimension if its not already one (like in the case of selecting a dataset via isel)
    This also allows you to loop through the time dimension, using the same code for a DataArray with one time step and a DataArray with many time steps 

    Args: 
        da (xr.DataArray)
        
    Returns: 
        da (xr.DataArray): da with time added as a coordinate and dimension, if it wasn't already 
    
    """
    
    if "time" not in da.dims: 
        if "time" not in da.coords: # If time isn't even a coordinate for some reason (like when taking the mean over the time dimension), set it as a coordinate here 
            da = da.assign_coords({"time":"unknown"}) 
        da = da.expand_dims("time")
    return da


def get_lat_lon_names(xr_ds): 
    """ Get name of lat and lon coordinate/dimension in xr_ds 
    I needed to do this because sometimes lat and lon is a coordinate while the horizontal dimensions is x,y. 
    The pcolormesh function will try to plot using x,y instead of lat,lon, which results in an empty plot 
    This automates the plotting function so the user doesn't have to inspect the dataframe/dataarray before plotting
    
    Args: 
        xr_ds (xr.Dataset or xr.DataArray)
    
    Returns: tuple
        lat_coord_name (str): name of latitude coordinate
        lon_coord_name (str): name of longitude coordinate
    
    """
    
    def get_coord_name(xr_ds, valid_names): 
        """ Check if xr.Dataset or xr.DataArray contains a coordinate or dimension value with the same name as a value in list valid_names 
        Used to check for name of lat and lon coordinate in a dataset. 
        For example, input valid_names=["latitude","longitudes","lat"] to see if xr_ds has a dim or coord value with the same name. The function will return "lat" if xr_ds contains it as a coord or dim

        Args: 
            xr_ds (xr.Dataset or xr.DataArray)

        Returns: 
            coord_str (str): whatever str value in valid_names list that is contained in xr_ds. If more than one value exists, the first will be returned

        """

        name_in_dims = [i for i in valid_names if i in list(xr_ds.dims) + list(xr_ds.coords)]
        if len(name_in_dims) == 0: 
            raise ValueError("Input dataset has no valid coordinate or dimension for latitude and/or longitude which is neccessary for plotting. \
                              \nInput dataset dimensions: "+", ".join(list(xr_ds.dims))+ \
                              "\nInput dataset coordinates: "+", ".join(list(xr_ds.coords))+ \
                               "\nAcceptable values: "+", ".join(valid_names))
        else: 
            coord_str = name_in_dims[0]
            return coord_str
    
    # Valid names for lat and lon coordinates 
    lat_list = ["latitude","longitudes","lat","lats","Latitude","Latitudes","Lat","Lats","LATITUDE","LATITUDES","LAT","LATS"]
    lon_list = ["longitude","longitudes","lon","lons","Longitude","Longitudes","Lon","Lons","LONGITUDE","Longitudes","LON","LONS"]
    
    # Get names of coordinates 
    lat_coord_name = get_coord_name(xr_ds, valid_names=lat_list)
    lon_coord_name = get_coord_name(xr_ds, valid_names=lon_list)
    
    # Return as a tuple 
    return lat_coord_name, lon_coord_name


def polarComparisonMaps(data1, data2, label1="data1", label2="data2", cmap="viridis", vmin=None, vmax=None, vmin_diff=None, vmax_diff=None, hemisphere="nh", lat=50, force_shared_timesteps=False, title=None): 
    """ Generate comparison maps with data1 on the left, data2 in the middle, and gridcell difference (data2 - data2) on the right. 
    This is a beefy function. It's set up to handle (almost) everything you throw at it, which is why it has a lot of complexity. 
    Older iterations of this function were a lot simpler, but they failed when I didn't give it the exact right thing for data1 and data2, so I added more complexity to this version.
    
    Args:
        data1 (xr.DataArray or xr.Dataset): data to plot on leftmost axis
        data2 (xr.DataArray or xr.Dataset): data to plot on middle axis
        label1 (str, optional): label to give axis title for data1 plot (default to "data1")
        label2 (str, optional): label to give axis title for data2 plot (default to "data2")
        cmap (str, optional): colormap to use for colorbar (default to "viridis")
        vmin (float, optional): minimum value to use for plotting data range (default to 1st percentile of either data1 or data2, whichever is smaller)
        vmax (float, optional): maximum value to use for plotting data range for model/obs data (default to 99th percentile of either data1 or data2, whichever is larger)
        vmin_diff (float, optional): minimum value to use for plotting data range (default to mapping default, with 0 at center)
        vmax_diff (float, optional): maximum value to use for plotting data range (default to mapping default, with 0 at center)
        hemisphere ('nh' or 'sh', optional): string indicating projection to use for map (default to 'nh')
        lat (float, optional): positive float value to restrict above (northern hemisphere) or below (southern hemisphere) (default to 50)
        force_shared_timesteps (bool, optional): Ensure that plots are only generated for shared timesteps? For example, set to False if you want to show the difference betweens two different months. (default to False)
        title (str, optional): main title to give figure (default to None)
        
    Returns: 
        plots displayed in notebook
    
    """
    
    # Add time as dimension and coordinate, if it's not already. See function add_time_dim for more info on why this is neccessary. Part of beefing up function.
    data1 = add_time_dim(data1)
    data2 = add_time_dim(data2)
    
    # Check if datasets have timesteps in common 
    if force_shared_timesteps == True: 
        shared_timesteps = [timestep for timestep in data1.time.values if timestep in data2.time.values]
        if len(shared_timesteps) == 0: # If no common timesteps are found, print warning and exit function 
            print("Input datasets do not contain any shared time values for which to plot. If you want to plot the data anyways, set force_shared_timesteps=False")
            print("Ending function. Returning None")
            return None
        else: 
            data1 = data1.sel(time=shared_timesteps)
            data2 = data2.sel(time=shared_timesteps)
    
    # If the user doesn't give a label, use the dataset's name instead of the default "data1", which isn't very descriptive
    if data1.name is not None and label1 == "data1": 
        label1 = data1.name
    if data2.name is not None and label2 == "data2": 
        label2 = data2.name
        
    # Get name of horizontal (lat, lon) coordinates. See function get_lat_lon_names for more information on why this was added  
    lat_coord_name1, lon_coord_name1 = get_lat_lon_names(data1)
    lat_coord_name2, lon_coord_name2 = get_lat_lon_names(data2)
    
    # Raise warning if lat and lon coord names are different for data1 and data2
    if (lat_coord_name1 != lat_coord_name2) or (lon_coord_name1 != lon_coord_name2): 
        print("WARNING: Horizontal coordinates are different for data1 and data2. You may see a blank plot for the differences unless you change the coordinate names.")

    # Compute min and max for plotting
    vmin_data1, vmax_data1 = compute_vmin_vmax(data1)
    vmin_data2, vmax_data2 = compute_vmin_vmax(data2) 
    vmin = vmin if vmin is not None else round(min(vmin_data1, vmin_data2),1) # Set to smallest value of the two 
    vmax = vmax if vmax is not None else round(max(vmax_data1, vmax_data2),1) # Set to largest value of the two 

    # Loop through each time coordinate, and plot! 
    for time1, time2 in zip(data1.time.values, data2.time.values): 

        # Select timestep 
        d1_timei = data1.sel(time=time1)
        d2_timei = data2.sel(time=time2)
        
        # Set up figure and axes 
        fig, axes = polar_fig_and_axes(hemisphere=hemisphere, figsize=(16,16), lat=lat)
        ax1, ax2, ax3 = axes[0], axes[1], axes[2]

        # Plot data for data1 and data 2
        im1 = d1_timei.plot.pcolormesh(ax=ax1, x=lon_coord_name1, y=lat_coord_name1, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=8, cmap=cmap, cbar_kwargs={'orientation':'horizontal', 'pad':0.02, 'extend':'both'})
        im2 = d2_timei.plot.pcolormesh(ax=ax2, x=lon_coord_name2, y=lat_coord_name2, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=8, cmap=cmap, cbar_kwargs={'orientation':'horizontal', 'pad':0.02, 'extend':'both'})

        # Compute and plot gridcell difference on 3rd axis 
        gridcell_diff = (d2_timei-d1_timei)
        im3 = gridcell_diff.plot.pcolormesh(ax=ax3, x=lon_coord_name1, y=lat_coord_name1, transform=ccrs.PlateCarree(), vmin=vmin_diff, vmax=vmax_diff, zorder=8, cmap='coolwarm', cbar_kwargs={'orientation':'horizontal', 'pad':0.02, 'extend':'both', 'label':'gridcell difference'})

        # Add titles to each axis 
        ax1.set_title(label1, fontdict={"fontsize":13})
        ax2.set_title(label2, fontdict={"fontsize":13})
        ax3.set_title("Griddcell difference ("+label2+" - "+label1+")", fontdict={"fontsize":13})
        
        # Add suptitle 
        if title is not None: 
            fig.suptitle(title, y=0.53, rotation='horizontal', fontsize=18)
        
        elif force_shared_timesteps == True: # Add timestep as title
            try: 
                time_str = "time = "+pd.to_datetime(time1).strftime("%Y-%m-%d") if time1 != "unknown" else "time = "+time1 # Time label for fig title
                fig.suptitle(time_str, y=0.53, rotation='horizontal', fontsize=18)
            except: 
                pass
        
        # Adjust subplots and display plot 
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.03, hspace=None) 
        plt.tight_layout()
        plt.show()


# -
def interactive_lineplot(df, title=None, grid=True, ylabel=None, xlabel=None, width=700, height=350, legend="right"): 
    """ Generate an interactive bokeh line plot
    
    Args: 
        df (pd.DataFrame): data to plot 
        grid (bool, optional): add gridlines to plot? (default to True)
        ylabel (str, optional): label for yaxis (default to "value")
        xlabel (str, optional): label for xaxis (default to "index")
        legend (str, optional): location for legend "top","bottom","right","left" (default to "right")
        title (str, optional): title to give plot (default to None)
        width (int, optional): width of plot (default to 650)
        height (int, optional): height of plot (default to 350)
    
    Returns: 
        lineplot (holoviews.core.overlay.Overlay)
    
    """
    
    lineplot = df.hvplot.line(grid=grid, legend=legend, ylabel=ylabel, xlabel=xlabel, title=title, width=width, height=height) \
                              * df.hvplot.scatter(marker='o') # Overlay scatter plot to add markers 
    return lineplot


def interactive_map(da, y="latitude", x="longitude", colorbar=True, vmin=None, vmax=None, clabel="", cmap="viridis", title=None, hemisphere="nh", lat=60, width=650, height=350): 
    """ Generate an interactive bokeh map. If the input da has more than one time value, a fun slider will be generated! 
    
    Args: 
        da (xr.DataArray): data to plot 
        y (str, optional): name of latitude dimension (default to "latitude")
        x (str, optional): name of longitude dimension (default to "longitude")
        colorbar (bool, optional): include colorbar? (default to True)
        vmin (float, optional): minumum for colorbar limit (default to 1st percentile of da. See compute_vmin_vmax function)
        vmax (float, optional): maximum for colorbar limit (default to 99th percentile da. See compute_vmin_vmax function)
        clabel (str, optional): label for colorbar (default to empty string)
        cmap (str, optional): colormap to use for data 
        title (str, optional): title to give plot (default to time dimension)
        hemisphere ('nh' or 'sh', optional): string indicating projection to use for map (default to 'nh')
        lat (float, optional): positive float value to restrict above (northern hemisphere) or below (southern hemisphere) (default to 60)
        width (int, optional): width of plot (default to 650)
        height (int, optional): height of plot (default to 350)
    
    Returns: 
        interactive_map (holoviews.core.spaces.DynamicMap)
    
    """
    
    if hemisphere=="nh": 
        max_lat = 90
        projection = ccrs.NorthPolarStereo(central_longitude=-45)
    elif hemisphere=="sh": 
        max_lat = -90 
        projection = ccrs.SouthPolarStereo()
        if lat > 0: 
            lat = -lat
    
    vmin_auto, vmax_auto = compute_vmin_vmax(da)
    clim = [vmin,vmax]
    if vmin is None: 
        clim[0] = vmin_auto
    if vmax is None: 
        clim[1] = vmax_auto 
    
    interactive_map = da.hvplot.quadmesh(y=y, x=x, 
                                         colorbar=colorbar, clim=tuple(clim), clabel=clabel, cmap=cmap, 
                                         title=title, features=["coastline"],
                                         projection=projection, geo=True, project=True, 
                                         width=width, height=height, ylim=(lat,max_lat))
    return interactive_map
