# +
""" plotting_utils.py 

Helper functions for generating maps and plots 

"""

import xarray as xr
import numpy as np 
import numpy.ma as ma
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from textwrap import wrap
import hvplot.xarray
import holoviews as hv
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh # Helps avoid some weird issues with the polar projection 


# -

def get_winter_data(da, year_start=None, start_month="Sep", end_month="Apr", force_complete_season=False):
    """ Select data for winter seasons corresponding to the input time range 
    
    Args: 
        da (xr.Dataset or xr.DataArray): data to restrict by time; must contain "time" as a coordinate 
        year_start (str, optional): year to start time range; if you want Sep 2019 - Apr 2020, set year="2019" (default to the first year in the dataset)
        start_month (str, optional): first month in winter (default to September)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
        
    Returns: 
        da_winter (xr.Dataset or xr.DataArray): da restricted to winter seasons 
    
    """
    if year_start is None: 
        print("No start year specified. Getting winter data for first year in the dataset")
        year_start = str(pd.to_datetime(da.time.values[0]).year)
    
    start_timestep = start_month+" "+str(year_start) # mon year 
    end_timestep = end_month+" "+str(int(year_start)+1) # mon year
    winter = pd.date_range(start=start_timestep, end=end_timestep, freq="MS") # pandas date range defining winter season
    months_in_da = [mon for mon in winter if mon in da.time.values] # Just grab months if they correspond to a time coordinate in da

    if len(months_in_da) > 0: 
        if (force_complete_season == True) and (all([mon in da.time.values for mon in winter])==False): 
            da_winter = None
        else: 
            da_winter = da.sel(time=months_in_da)
    else: 
        da_winter = None
        
    return da_winter


def compute_gridcell_winter_means(da, years=None, start_month="Nov", end_month="Apr", force_complete_season=False): 
    """ Compute winter means over the time dimension. Useful for plotting as the grid is maintained. 
    
    Args: 
        da (xr.Dataset or xr.DataArray): data to restrict by time; must contain "time" as a coordinate 
        years (list of str): years over which to compute mean (default to unique years in the dataset)
        year_start (str, optional): year to start time range; if you want Nov 2019 - Apr 2020, set year="2019" (default to the first year in the dataset)
        start_month (str, optional): first month in winter (default to November)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
    
    Returns: 
        merged (xr.DataArray): DataArray with winter means as a time coordinate
    """
    
    if years is None: 
        years = np.unique(pd.to_datetime(da.time.values).strftime("%Y")) # Unique years in the dataset 

    winter_means = []
    for year in years: # Loop through each year and grab the winter months, compute winter mean, and append to list 
        da_winter_i = get_winter_data(da, year_start=year, start_month=start_month, end_month=end_month, force_complete_season=force_complete_season)
        if da_winter_i is None: 
            continue
        da_mean_i = da_winter_i.mean(dim="time", keep_attrs=True) # Comput mean over time dimension

        # Assign time coordinate 
        time_arr = pd.to_datetime(da_winter_i.time.values)
        da_mean_i = da_mean_i.assign_coords({"time":time_arr[0].strftime("%b %Y")+" - "+time_arr[-1].strftime("%b %Y")})
        da_mean_i = da_mean_i.expand_dims("time")

        winter_means.append(da_mean_i)

    merged = xr.merge(winter_means) # Combine each winter mean Dataset into a single Dataset, with the time period maintained as a coordinate
    merged = merged[list(merged.data_vars)[0]] # Convert to DataArray
    merged.time.attrs["description"] = "Time period over which mean was computed" # Add descriptive attribute 
    return merged 


def staticArcticMaps(da, title=None, out_str="out", cmap="viridis", col=None, col_wrap=3, vmin=None, vmax=None, set_cbarlabel = '', min_lat=50, savefig=True): 
    """ Show data on a basemap of the Arctic. Can be one month or multiple months of data. 
    Creates an xarray facet grid. For more info, see: http://xarray.pydata.org/en/stable/user-guide/plotting.html
    
    Args: 
        da (xr DataArray): data to plot
        title (str, optional): title string for plot
        out_str (str, optional): output string when saving
        cmap (str, optional): colormap to use (default to viridis)
        col (str, optional): coordinate to use for creating facet plot (default to "time")
        col_wrap (int, optional): number of columns of plots to display (default to 3, or None if time dimension has only one value)
        vmin (float, optional): minimum on colorbar (default to 1st percentile)
        vmax (float, optional): maximum on colorbar (default to 99th percentile)
        min_lat (float, optional): minimum latitude to set extent of plot (default to 50 deg lat)
        set_cbarlabel (str, optional): set colorbar label
        savefig (bool): output figure
    
    Returns:
        Figure displayed in notebook 
    
    """ 
    # Compute min and max for plotting
    def compute_vmin_vmax(da): 
        vmin = np.nanpercentile(da.values, 1)
        vmax = np.nanpercentile(da.values, 99)
        return vmin, vmax
    vmin_data, vmax_data = compute_vmin_vmax(da)
    vmin = vmin if vmin is not None else vmin_data # Set to smallest value of the two 
    vmax = vmax if vmax is not None else vmax_data # Set to largest value of the two 
    
    # All of this col and col_wrap maddness is to try and make this function as generalizable as possible
    # This allows the function to work for DataArrays with multiple coordinates, different coordinates besides time, etc! 
    if col is None: 
        col = "time"
        try: # Assign time coordinate if it doesn't exist
            da["time"]
        except AttributeError: 
            da = da.assign_coords({col:"unknown"})
    col = col if sum(da[col].shape) > 1 else None
    if col is not None: 
        if sum(da[col].shape)<=1: 
            col_wrap = None
    
    # Plot
    if len(set_cbarlabel)==0:
        set_cbarlabel=da.attrs["long_name"]+' ['+da.attrs["units"]+']'

    im = da.plot(x="longitude", y="latitude", col_wrap=col_wrap, col=col, transform=ccrs.PlateCarree(), cmap=cmap, zorder=8, 
             cbar_kwargs={'pad':0.02,'shrink': 0.8,'extend':'both', 'label':set_cbarlabel},
             vmin=vmin, vmax=vmax, 
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
        ax.set_extent([-179, 179, min_lat, 90], crs=ccrs.PlateCarree()) # Set extent to zoom in on Arctic
    
    # Get figure
    fig = plt.gcf()
    
    # Set title 
    if (sum(ax_iter.shape) == 0) and (title is not None): 
        ax.set_title(title, fontsize=12, horizontalalignment="center", x=0.45, y=1.06, fontweight='medium')
    elif title is not None:
        fig.suptitle(title, fontsize=12, horizontalalignment="center", x=0.45, y=1.06, fontweight='medium')
    
    # save figure
    if savefig:
        plt.savefig('./figs/maps_'+out_str+'.png', dpi=400, facecolor="white", bbox_inches='tight')

    plt.close() # Close so it doesnt automatically display in notebook 
    return fig


def staticArcticMaps_overlayDrifts(da, drifts_x, drifts_y, alpha=1, vector_val=0.1, scale_vec=0.5, res=6, units_vec=r'm s$^{-1}$', title=None, out_str="out", cmap="viridis", col=None, col_wrap=3, vmin=None, vmax=None, set_cbarlabel = '', min_lat=50, savefig=True, figsize=(6,6)): 
    """ Show data on a basemap of the Arctic. Can be one month or multiple months of data. Overlay drift vectors on top 
    Creates an xarray facet grid. For more info, see: http://xarray.pydata.org/en/stable/user-guide/plotting.html
    
    Args: 
        da (xr DataArray): data to plot
        drifts_x (xr.DataArray): sea ice drifts along-x component of the ice motion
        drifts_y (xr.DataArray): sea ice drifts along-y component of the ice motion
        alpha (float 0-1, optional): Set this variable if you want da to have a reduced opacity (default to 1)
        res (int, optional): resolution of vectors (default to 6; plot 1 out of every 6 vectors)
        title (str, optional): title string for plot
        out_str (str, optional): output string when saving
        cmap (str, optional): colormap to use (default to viridis)
        col (str, optional): coordinate to use for creating facet plot (default to "time")
        col_wrap (int, optional): number of columns of plots to display (default to 3, or None if time dimension has only one value)
        vmin (float, optional): minimum on colorbar (default to 1st percentile)
        vmax (float, optional): maximum on colorbar (default to 99th percentile)
        min_lat (float, optional): minimum latitude to set extent of plot (default to 50 deg lat)
        set_cbarlabel (str, optional): set colorbar label
        savefig (bool): output figure
    
    Returns:
        Figure displayed in notebook 
    
    """ 
    # Make sure alpha is between 0 and 1 
    if alpha > 1: 
        print("Argument alpha must be between 0 and 1. You inputted " +str(alpha)+ ". Setting alpha to 1.")
        alpha = 1 
    elif alpha < 0: 
        print("Argument alpha must be between 0 and 1. You inputted " +str(alpha)+ ". Setting alpha to 0.5.")
        alpha = 0.5
    elif alpha == 0: 
        print("You set alpha=0. This indicates full transparency of the input data. No data will be displayed on the map.")
    
    # Check that drifts and da have the same time coordinates 
    for drift in [drifts_x,drifts_y]:
        equality = (da.time.values  == drift.time.values)
        if type(equality) == np.ndarray:
            if not all(equality): 
                raise ValueError("Drifts vectors and input DataArray must have the same time coordinates")
        elif (equality==False):
            raise ValueError("Drifts vectors and input DataArray must have the same time coordinates")

    # Compute min and max for plotting
    def compute_vmin_vmax(da): 
        vmin = np.nanpercentile(da.values, 1)
        vmax = np.nanpercentile(da.values, 99)
        return vmin, vmax
    vmin_data, vmax_data = compute_vmin_vmax(da)
    vmin = vmin if vmin is not None else vmin_data # Set to smallest value of the two 
    vmax = vmax if vmax is not None else vmax_data # Set to largest value of the two 
    
    # All of this col and col_wrap maddness is to try and make this function as generalizable as possible
    # This allows the function to work for DataArrays with multiple coordinates, different coordinates besides time, etc! 
    if col is None: 
        col = "time"
        try: # Assign time coordinate if it doesn't exist
            da["time"]
        except AttributeError: 
            da = da.assign_coords({col:"unknown"})
    col = col if sum(da[col].shape) > 1 else None
    if col is not None: 
        if sum(da[col].shape)<=1: 
            col_wrap = None
            
    # Plot
    if len(set_cbarlabel)==0:
        set_cbarlabel=da.attrs["long_name"]+' ['+da.attrs["units"]+']'

    im = da.plot(x="longitude", y="latitude", col_wrap=col_wrap, col=col, transform=ccrs.PlateCarree(), cmap=cmap, 
                 cbar_kwargs={'pad':0.02,'shrink': 0.8,'extend':'both', 'label':set_cbarlabel},
                 vmin=vmin, vmax=vmax, zorder=2, alpha=alpha, 
                 subplot_kws={'projection':ccrs.NorthPolarStereo(central_longitude=-45)})
    
    # Iterate through axes and add features 
    ax_iter = im.axes
    if type(ax_iter) != np.array: # If the data is just a single month, ax.iter returns an axis object. We need to iterate through a list or array
        ax_iter = np.array(ax_iter)
    
    i = 0
    try: 
        num_timesteps = len(da.time.values)
    except: 
        num_timesteps = 1
    for ax, i in zip(ax_iter.flatten(), range(num_timesteps)):

            # Add drifts 
            if num_timesteps == 1: 
                drifts_xi = drifts_x.copy()
                drifts_yi = drifts_y.copy()
            else: 
                drifts_xi = drifts_x.isel(time=i).copy()
                drifts_yi = drifts_y.isel(time=i).copy()
            Q = ax.quiver(drifts_x.xgrid[::res, ::res], drifts_y.ygrid[::res, ::res], 
                          ma.masked_where(np.isnan(drifts_xi[::res, ::res]), drifts_xi[::res, ::res]),
                          ma.masked_where(np.isnan(drifts_yi[::res, ::res]), drifts_yi[::res, ::res]) , units='inches', scale=scale_vec, zorder=10)
            ax.quiverkey(Q, 0.85, 0.88, vector_val, str(vector_val)+' '+units_vec, coordinates='axes', zorder=11)   

            ax.coastlines(linewidth=0.15, color = 'black', zorder = 8) # Coastlines
            ax.add_feature(cfeature.LAND, color ='0.95', zorder = 5)    # Land
            ax.add_feature(cfeature.LAKES, color = 'grey', zorder = 5)  # Lakes
            ax.gridlines(draw_labels=False, linewidth=0.25, color='gray', alpha=0.7, linestyle='--', zorder=6) # Gridlines
            ax.set_extent([-179, 179, min_lat, 90], crs=ccrs.PlateCarree()) # Set extent to zoom in on Arctic
        
    # Get figure
    fig = plt.gcf()
    
    # Set title 
    if (sum(ax_iter.shape) == 0) and (title is not None): 
        ax.set_title(title, fontsize=12, horizontalalignment="center", x=0.45, y=1.06, fontweight='medium')
    elif title is not None:
        fig.suptitle(title, fontsize=12, horizontalalignment="center", x=0.45, y=1.06, fontweight='medium')
    
    # save figure
    if savefig:
        plt.savefig('./figs/maps_'+out_str+'.png', dpi=400, facecolor="white", bbox_inches='tight')
        
    plt.close() # Close so it doesnt automatically display in notebook 
    return fig


def interactiveArcticMaps(da, clabel=None, cmap="viridis", colorbar=True, vmin=None, vmax=None, title="", ylim=(60,90), frame_width=250, slider=True, cols=3): 
    """ Generative one or more interactive maps 
    Using the argument "slide", the user can set whether each map should be displayed together, or displayed in the form of a slider 
    To show each map together (no slider), set slider=False
    
    Args: 
        da (xr.Dataset or xr.DataArray): data 
        clabel (str, optional): colorbar label (default to "long_name" and "units" if given in attributes of da)
        cmap (str, optional): matplotlib colormap to use (default to "viridis")
        colorbar (bool, optional): show colorbar? (default to True)
        vmin (float, optional): minimum on colorbar (default to 1st percentile)
        vmax (float, optional): maximum on colorbar (default to 99th percentile)
        title (str, optional): main title to give plot (default to no title)
        ylim (tuple, optional): limits of yaxis in the form min latitude, max latitude (default to (60,90))
        frame_width (int, optional): width of frame. sets figure size of each map (default to 250)
        slider (bool, optional): if da has more than one time coordinate, display maps with a slider? (default to True)
        cols (int, optional): how many columns to show before wrapping, if da has more than one time coordinate (default to 3)
    
    Returns: 
        pl (Holoviews map)
    
    """
    # Compute min and max for plotting
    def compute_vmin_vmax(da): 
        vmin = np.nanpercentile(da.values, 1)
        vmax = np.nanpercentile(da.values, 99)
        return vmin, vmax
    vmin_data, vmax_data = compute_vmin_vmax(da)
    vmin = vmin if vmin is not None else vmin_data # Set to smallest value of the two 
    vmax = vmax if vmax is not None else vmax_data # Set to largest value of the two 
    
    #https://hvplot.holoviz.org/user_guide/Subplots.html
    subplots=False
    shared_axes=False
    show_title=True
    if ("time" in da.coords):
        if (sum(da["time"].shape) > 1): 
            subplots=True
            shared_axes=True
            if slider==True and title=="": 
                show_title=False # We don't want to remove the title for the slider plots since it removes the time from the title 
        
    if clabel is None and ("long_name" in da.attrs): # Add a logical colorbar label 
        clabel=da.attrs["long_name"]
        if "units" in da.attrs: 
            clabel+=" ("+da.attrs["units"]+")"
        
    pl = da.hvplot.quadmesh(y="latitude", x="longitude",
                            projection=ccrs.NorthPolarStereo(central_longitude=-45), 
                            features=["coastline"], # Add coastlines 
                            colorbar=colorbar, clim=(vmin,vmax), cmap=cmap, clabel=clabel, # Colorbar settings 
                            project=True, ylim=ylim, frame_width=frame_width,
                            subplots=subplots, shared_axes=shared_axes,
                            dynamic=False) 
    if slider==False: # Set number of columns 
        pl = pl.layout().cols(cols)
    
    if show_title==True: 
        pl.opts(title=title) # Add title
    hv.output(widget_location="bottom")
    return pl 


def interactive_winter_mean_maps(da, years=None, end_year=None, start_month="Sep", end_month="Apr", force_complete_season=False, clabel=None, cmap="viridis", colorbar=True, vmin=0, vmax=4, title="", ylim=(60,90), frame_width=250, slider=True, cols=3): 
    """ Generate interactive maps of winter mean data 
    Note: this function builds off the functions get_winter_data and interactiveArcticMaps.
    
    Args: 
        da (xr.Dataset or xr.DataArray): data; must contain "time" coordinate
        years (list of str): years over which to compute mean (default to unique years in the dataset)
        start_month (str, optional): first month in winter (default to September)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
        clabel (str, optional): colorbar label (default to "long_name" and "units" if given in attributes of da)
        cmap (str, optional): matplotlib colormap to use (default to "viridis")
        colorbar (bool, optional): show colorbar? (default to True)
        vmin (float, optional): minimum on colorbar (default to 0)
        vmax (float, optional): maximum on colorbar (default to 4)
        title (str, optional): main title to give plot (default to no title)
        ylim (tuple, optional): limits of yaxis in the form min latitude, max latitude (default to (60,90))
        frame_width (int, optional): width of frame. sets figure size of each map (default to 250)
        slider (bool, optional): if da has more than one time coordinate, display maps with a slider? (default to True)
        cols (int, optional): how many columns to show before wrapping, if da has more than one time coordinate (default to 3)
    
    Returns: 
        pl_means (Holoviews map)
    
    """
    
    winter_means_da = compute_gridcell_winter_means(da, years=years, start_month=start_month, end_month=end_month, force_complete_season=force_complete_season)

    pl_means = interactiveArcticMaps(winter_means_da, 
                                    clabel=clabel, cmap=cmap, colorbar=colorbar, 
                                    vmin=vmin, vmax=vmax, title=title, 
                                    ylim=ylim, frame_width=frame_width, slider=slider, cols=cols)
    hv.output(widget_location="bottom")
    return pl_means


def static_winter_comparison_lineplot(da, da_unc=None, years=None, figsize=(5,3), start_month="Sep", end_month="Apr", title="", set_ylabel = '', set_units = '', legend=True, savefig=True, save_label='', force_complete_season=False): 
    """ Make a lineplot with markers comparing monthly mean data across winter seasons 
    
    Args: 
        da (xr.DataArray): data to plot and compute mean for; must contain "time" as a coordinate 
        years (list of str): list of years for which to plot data. 2020 would correspond to the winter season defined by start month 2020 - end month 2021 (default to all unique years in da)
        title (str, optional): title to give plot (default to no title) 
        set_ylabel (str, optional): prescribed y label string
        set_units (str, optional): prescribed y label unit string
        legend (bool): print legend
        savefig (bool): output figure
        save_label (str, optional): additional string for output
        figsize (tuple, optional): figure size to display in notebook (default to (5,3))
        start_month (str, optional): first month in winter (default to September)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
        
       Returns: 
           Figure displayed in notebook
        
    """
    if years is None: 
        years = np.unique(pd.to_datetime(da.time.values).strftime("%Y")) # Unique years in the dataset 
        print("No years specified. Using "+", ".join(years))
    
    # Set up x-axis 
    # This avoids having a set x-axis of winter months between Sep-Apr, even if there's no data for Sep, Oct etc 
    yr = 2000 
    if end_month not in ["Oct","Nov","Dec"]: 
        yr_end = yr+1
    else: 
        yr_end = yr
    xaxis_months = pd.date_range(start_month+"-"+str(yr), end_month+"-"+str(yr_end), freq="M").strftime("%b")
    
    # Set up plot 
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(xaxis_months, np.empty((len(xaxis_months),1))*np.nan, color=None, label=None) # Set x axis using winter months 
    gridlines = plt.grid(b = True, linestyle = '-', alpha = 0.2) # Add gridlines 


    fmts = ['mo-','cs-','yv-','b*-.','r.-','gD--','k2-.']
    for year, fmt in zip(years, fmts*100): 
        winter_da = get_winter_data(da, year_start=year, start_month=start_month, end_month=end_month, force_complete_season=force_complete_season) # Get data from that winter 
        if winter_da is None: # In case the user inputs a year that doesn't have data, skip this loop iteration to avoid appending None
            continue
        y = winter_da.mean(dim=["x","y"], keep_attrs=True)
        x = pd.to_datetime(y.time.values)
        ax.plot(x.strftime("%b"), y, fmt, label="Winter "+str(x.year[0])+"-"+str(x.year[-1])[2:])

        if da_unc is not None:
            # Get uncertaintiy data from that winter 
            winter_da_unc = get_winter_data(da_unc, year_start=year, start_month=start_month, end_month=end_month, force_complete_season=force_complete_season) 
            if winter_da_unc is None: # In case the user inputs a year that doesn't have data, skip this loop iteration to avoid appending None
                continue
            yu = winter_da_unc.mean(dim=["x","y"], keep_attrs=True)    
            ax.fill_between(x.strftime("%b"), y - yu, y + yu, facecolor = fmt[0], alpha = 0.1, edgecolor = 'none')
    

    # Add legend, title, and axis labels, and display plot in notebook 
    if legend:
        plt.legend(fontsize=8, loc="best")
    
    plt.title(title, fontsize=9)
    if len(set_ylabel)>0:
        ylabel=set_ylabel
    elif "long_name" in da.attrs: 
        ylabel = da.attrs["long_name"]
        if "units" in da.attrs: 
            ylabel+=" ("+da.attrs["units"]+")"
        ylabel="\n".join(wrap(ylabel, 35))
    else: 
        ylabel=None

    plt.ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='both', which='major', labelsize=8)
   
   # reduce white space
    plt.tight_layout()

    # save figure
    if savefig:
        plt.savefig('./figs/'+da.attrs["long_name"]+start_month+end_month+str(years[0])+'-'+str(years[-1])+save_label+'.pdf', 
                    dpi=300, facecolor="white", bbox_inches='tight')

    plt.show()


def interactive_winter_comparison_lineplot(da, years=None, title="Winter comparison", frame_width=600, frame_height=350, start_month="Sep", end_month="Apr", force_complete_season=False):
    """ Make a bokeh lineplot with markers comparing monthly mean data across winter seasons 
    
    Args: 
        da (xr.DataArray): data; must contain "time" coordinate
        years (list of str): list of years for which to plot data. 2020 would correspond to the winter season defined by start month 2020 - end month 2021 (default to all unique years in da)
        title (str, optional): title to give plot (default to "Winter comparison") 
        frame_width (int, optional): width of figure (default to 600) 
        frame_height (int, optional): height of figure (default to 350) 
        start_month (str, optional): first month in winter (default to September)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
        
       Returns: 
           pl (bokeh lineplot) 
        
    """
    
    if years is None: 
        years = np.unique(pd.to_datetime(da.time.values).strftime("%Y")) # Unique years in the dataset 
    
    winter_means_list = []
    for year in years:
        winter_da = get_winter_data(da, year_start=year, start_month=start_month, end_month=end_month, force_complete_season=force_complete_season) # Get data from that winter 
        if winter_da is None: # In case the user inputs a year that doesn't have data, skip this loop iteration to avoid appending None
            continue
        winter_means_list.append(winter_da)
            
    # Sort by longest --> shortest. This avoids weird issues with x axis trying to be in time order 
    winter_means_list_sorted = sorted(winter_means_list, key=lambda l: (len(l), l))[::-1]

    # Combine plots and display
    i = 0
    for da_sorted in winter_means_list_sorted: 
        winter_mean = da_sorted.mean(dim=["x","y"], keep_attrs=True) # Compute mean 
        winter_mean["time"] = pd.to_datetime(da_sorted["time"].values).strftime("%b") # Reassign the time coordinate to be just the months (Nov, Dec, ect). This allows you to easily overlay the plots on top of each other, since they share an axis
        time_str = pd.to_datetime(da_sorted.time).strftime("%Y") # Get time coordinate as string value

        pl = winter_mean.hvplot(grid=True, label="Winter "+time_str[0]+"-"+time_str[-1], frame_width=frame_width, frame_height=frame_height) * winter_mean.hvplot.scatter(marker='o') # Overlay scatter plot to add markers
        if i == 0:
            pl_tot = pl
        else: 
            pl_tot *= pl 
        i+=1
        
    winters_all = pl_tot.opts(hv.opts.Layout(shared_axes=True, merge_tools=True)) # Combine lineplots into a single figure 
    winters_all.opts(title=title) # Add a title 
    return winters_all
