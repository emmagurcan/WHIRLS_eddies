#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 11:24:24 2025

@author: xmiram
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy.signal import fftconvolve
import gsw 
import datetime as dt

import glidertools as gt

from matplotlib import rc
rc('font', weight='normal',size=16)

# Path to save figures -> insert your path here
path_fig = '/Users/xmiram/Documents/Postdoc/Gough/Figures/'

#%%
# =============================================================================
#        Functions
# =============================================================================

# Get dates from unix time
def get_date(self):
    return [dt.datetime.utcfromtimestamp(i) for i in self]

# My favorite smoothing function
def smoo(self,nw):
	win = np.blackman(nw)
	#win=win[None,:]*win[:,None] #for 2D plots
	peso= fftconvolve(np.ones(self.shape),win,mode='same')
	dataout = fftconvolve(self,win,mode='same')/peso
	return dataout

# Function to plot all lines and collect handles
def plot_profiles(current_ax, x_data, y_data, gl_x_col, is_main_plot=False):
    
    # Ship CTD (Black line) - Only labeled on the main plot for legend collection
    current_ax.plot(x_data, y_data, 'k', lw=4)
    ship_line, = current_ax.plot(x_data, y_data, '0.5', lw=3)
    
    # Glider Profiles (Colored lines)
    glider_handles = []
    for i, ind in enumerate(plot_dives):
        # Use a label only if it's the main plot for legend handles
        label_text = legend_times[i] if is_main_plot else None 
        
        glider_line, = current_ax.plot(
            gl[gl_x_col].where(gl.dives == ind), 
            -gl.ctd_pressure.where(gl.dives == ind), 
            c=dc[i], 
            lw=2, 
            label=label_text
        )
        if is_main_plot:
            glider_handles.append(glider_line)
    
    if is_main_plot:
        # Return ship handle and all glider handles
        return [ship_line] + glider_handles
    return None # Return None for inset/salinity plots

def closeto(x,number):
    ind = np.argwhere(np.abs(x-number)==np.nanmin(np.abs(x-number)))
    return ind[0][0]

def smart_extract_and_rename(df, keywords_with_std_names, fuzzy=True):
    """
    Extract columns matching keywords, rename to standard names, and store metadata.
    
    Parameters
    ----------
    df : pandas.DataFrame
    keywords_with_std_names : dict
        Mapping of keywords to standard variable names, e.g.,
        {"temperature": "temperature", "salinity": "salinity", "oxygen": "oxygen"}
    fuzzy : bool
        If True, allows partial matching.
    
    Returns
    -------
    xarray.Dataset with renamed variables and metadata
    """
    cols = df.columns.tolist()
    all_matches = {}
    
    for keyword, std_name in keywords_with_std_names.items():
        if fuzzy:
            matches = [c for c in cols if keyword.lower() in c.lower()]
        else:
            matches = [c for c in cols if c.lower() == keyword.lower()]
        
        if not matches:
            print(f"Warning: No variables found matching '{keyword}'")
            continue
        
        # Take the first match as representative (could be extended to multiple)
        matched_col = matches[0]
        all_matches[std_name] = matched_col
        print(f"Matched '{keyword}' -> '{matched_col}' as '{std_name}'")
    
    if not all_matches:
        raise ValueError("No variables matched any of the keywords.")
    
    # Extract columns
    ds = df[list(all_matches.values())].to_xarray()
    
    # Rename to standard names and store original name as metadata
    rename_dict = {}
    for std_name, original_name in all_matches.items():
        rename_dict[original_name] = std_name
        ds[original_name].attrs["original_name"] = original_name
        # Optionally, if units column exists, you can store units here too
        # ds[original_name].attrs["units"] = ...
    
    ds = ds.rename(rename_dict)
    return ds

def grid_seaglider(gl, bin_edges, exclude_vars=None):
    """
    Bin Seaglider data in pressure using the existing 'dives' profile labels.
    
    Parameters
    ----------
    gl : xarray.Dataset
        Must contain 'ctd_pressure' and 'dives'.
    bin_edges : array-like
        Pressure bin edges (e.g., np.arange(0, 1001, 1)).
    exclude_vars : list
        Variables that should not be binned.
    
    Returns
    -------
    xr.Dataset
        Binned data with dimensions: pressure x profile
    """

    import numpy as np
    import xarray as xr

    # -------------------------------------------------------
    # Basic checks
    # -------------------------------------------------------
    if "ctd_pressure" not in gl:
        raise KeyError("Dataset must contain 'ctd_pressure'.")
    if "dives" not in gl:
        raise KeyError("Dataset must contain 'dives'.")

    # Default excluded variables
    if exclude_vars is None:
        exclude_vars = ["ctd_time", "ctd_time_dt64"] #"longitude", "latitude"

    # Identify data variables to bin
    data_vars = [v for v in gl.data_vars if v not in exclude_vars and v not in ["ctd_pressure", "dives", "ctd_time_dt64"]]
    data_vars = data_vars + [v for v in gl.coords if v not in ["ctd_time_dt64"]]
    
    # Compute bin centers
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
    n_bins = len(bin_centers)

    # Unique dive/climb profile IDs
    profiles = np.unique(gl["dives"].values)
    n_profiles = len(profiles)

    # Convert to dataframe for easy grouping
    df = gl.to_dataframe().reset_index()

    # Pre-allocate output arrays
    binned = {
        var: np.full((n_bins, n_profiles), np.nan)
        for var in data_vars
    }

    # -------------------------------------------------------
    # Loop over each profile
    # -------------------------------------------------------
    for i, prof in enumerate(profiles):

        group = df[df["dives"] == prof]
        if group.empty:
            continue

        p = group["ctd_pressure"].values
        bin_idx = np.digitize(p, bin_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < n_bins)

        # Loop through pressure bins in this profile
        for b in np.unique(bin_idx[valid]):
            mask = bin_idx == b
            if np.any(mask):
                for var in data_vars:
                    binned[var][b, i] = np.nanmedian(group[var].values[mask])

    # -------------------------------------------------------
    # Build xarray output dataset
    # -------------------------------------------------------
    ds_out = xr.Dataset(
        data_vars={f"{var}_binned": (["pressure", "profile"], binned[var]) for var in data_vars},
        coords={
            "pressure": bin_centers,
            "profile": profiles,
        },
        attrs={
            "description": "Seaglider pressure-binned profiles",
            "binning": "median within each pressure bin and dive/climb profile"
        }
    )

    return ds_out

#%%
# =============================================================================
# Load SG data
# =============================================================================

# Testing with one file
# file = '/Users/xmiram/SG675_Gough_backup/p6750002.nc'
# test = xr.open_dataset(file)

# Grabbing all glider files
filenames = '/Users/xmiram/SG675_Gough_backup/p*.nc'

names = [
    'ctd_depth',
    'ctd_time',
    'ctd_pressure',
    'salinity',
    'temperature',
    'conductivity',
    #'wlbb2fl_sig470nm_adjusted',
    #'wlbb2fl_sig700nm_adjusted',
    'wlbb2fl_sig695nm_adjusted', # chlorophyll-a concentration, microg/L
    #'aanderaa4831_results_time',
    'aanderaa4831_instrument_dissolved_oxygen',
    #'aanderaa4831_drift_gain',
    #'aanderaa4831_dissolved_oxygen',
    #'aanderaa4831_dissolved_oxygen_qc',
    #'aanderaa4831_qc',
    'eng_aa4831_O2',
    #'eng_wlbb2fl_ref700nm',
    #'eng_wlbb2fl_sig700nm',
    #'eng_wlbb2fl_ref470nm',
    #'eng_wlbb2fl_sig470nm',
    #'eng_wlbb2fl_ref695nm',
    #'eng_wlbb2fl_sig695nm',
    'eng_qsp_PARuV', # PAR
    'gsw_sigma0'    
]

tt = gt.load.seaglider_basestation_netCDFs(
    filenames, names,
    return_merged=False,
    keep_global_attrs=False
)

gl = tt["sg_data_point"]

# Alternatively, you can load the entire time series...
# file = '/Users/xmiram/SG675_Gough_backup/sg675_SG675_Gough_SAMBA_timeseries.nc'
# gl = xr.open_dataset(file)

#%%
# =============================================================================
# Load CTD variables
# =============================================================================

# Load Excel file into a pandas DataFrame
file = '/Users/xmiram/SG_taskforce/sea_glider_processing/CTD_data/Gough2025_Glider_Station_CTD_profile_data.csv' 

df = pd.read_csv(file)

print(df.columns.tolist())

keywords = {
    "time": "time",
    "temperature": "T",
    "salinity": "S",
    "oxygen": "Oxy",
    "conductivity": "C",
    "pressure": "press",
    "fluorescence": "flu",
    "turbidity": "turb"
    
}

ctd = smart_extract_and_rename(df, keywords)

#%%

# =============================================================================
# TS profiles comparison, CTD and glider
# =============================================================================

# Choosing some colors for ploting
colormap_name = 'Reds'
num_colors = 6
cmap = plt.cm.get_cmap(colormap_name)
color_indices = np.linspace(0.3, 1, num_colors)
dc = [cmap(i) for i in color_indices]

# Define the dive indices to iterate over -> 13 is the first profile we have to 1000 m
plot_dives = [2, 2.5, 5, 8, 10, 13]

# ----------------------------------------
# Getting the lables for each profile -> time legend
legend_times = []
for ind in plot_dives:
    # 1. Get the minimum time as an xarray DataArray
    first_time_da = gl.where(gl.dives == ind).dropna("sg_data_point").ctd_time_dt64.mean()
    
    # 2. Extract the underlying NumPy scalar value
    time_scalar = first_time_da.values.item() 
    
    # 3. Check if the value is a valid time (i.e., not Not a Time/NaN)
    # np.isnat checks if a numpy datetime object is NaT
    if pd.notna(time_scalar):
        # Convert the numpy.datetime64/pandas NaT back to a pandas Timestamp for reliable formatting
        first_time_timestamp = pd.to_datetime(time_scalar)
        
        # Format the time (e.g., 'Sep 17, 10:00')
        legend_times.append("SG, "+first_time_timestamp.strftime('%b %d, %H:%M'))
    else:
        # Handle cases where the dive has no associated time (e.g., all times were null)
        legend_times.append('No Data Time')

# Combine ship label and time labels for the final legend
all_labels = ['Ship CTD, Sep 17 11:16'] + legend_times # Not sure if the date and time for the CTD is correct...
# ----------------------------------------

# Figure begins here! ----------------------------------------
fig, [ax, bx] = plt.subplots(1, 2, figsize=(10, 15), constrained_layout=True, sharey=True)

# --- ax: Temperature Profile (Main Plot) ---
# Pass is_main_plot=True to collect handles
all_handles = plot_profiles(ax, ctd.T, -ctd.press, 'temperature', is_main_plot=True)

ax.set_ylim(-1000,0) 
ax.set_xlim(3.5,17.2)
ax.set_yticks(np.arange(-1000,100,100),labels=np.arange(1000,-100,-100))
ax.set_xlabel(r'Temperature ($^{\circ}$C)')
ax.set_ylabel('Pressure (dbar)')
ax.grid(lw=0.2)

# Adding legend to ax (Temperature)
ax.legend(
    all_handles, 
    all_labels, 
    #title='Profile Start Time',
    loc='lower right', 
    bbox_to_anchor=(0.98, 0.02), # Near the bottom-right corner of the plot area
    fontsize=14,
    frameon=False,
    markerfirst=False
    
)

# --- ax inset: Temperature Zoom (0-100m) ---
ax_zoom = ax.inset_axes([0.000, 0.7, 0.5, 0.3]) 
plot_profiles(ax_zoom, ctd.T, -ctd.press, 'temperature')

ax_zoom.set_ylim(-100, 0)
ax_zoom.set_xlim(16.1, 17.2)
ax_zoom.tick_params(labelsize=8)
ax_zoom.grid(lw=0.2)
ax_zoom.set_yticks([])

# Create a connection between plots
ax.indicate_inset_zoom(ax_zoom,edgecolor="k",lw=2)

# --- bx: Salinity Profile (Main Plot) ---
plot_profiles(bx, ctd.S, -ctd.press, 'salinity') # No need to collect handles here
bx.set_xlabel(r'Salinity (PSU)')
bx.grid(lw=0.2)

# --- bx inset: Salinity Zoom (0-100m) ---
bx_zoom = bx.inset_axes([0.000, 0.7, 0.5, 0.3]) 
plot_profiles(bx_zoom, ctd.S, -ctd.press, 'salinity')

bx.indicate_inset_zoom(bx_zoom,edgecolor="k",lw=2)

bx_zoom.set_ylim(-100, 0)
bx_zoom.set_xlim(35.44, 35.55)
bx_zoom.set_xticks(np.arange(35.45, 35.56,0.02))
bx_zoom.tick_params(labelsize=8)
bx_zoom.grid(lw=0.2)
bx_zoom.set_yticks([])

fig.savefig(path_fig+'CTD_cal_TSprofiles_Gough_Melktert_Sep2025.png',dpi=300)

#%%

# =============================================================================
# Define the biases for T and S 
# =============================================================================

# 1. Bin the SG data in pressure to match the CTD pressure

# Data binning -----------------------------
bin_edges = np.arange(-0.5, 1002, 1)

ds_bin = grid_seaglider(gl, bin_edges)

# Uncomment this to take a look at the dataset
# ds_bin

# =============================================================================
# Temperature
# =============================================================================

# This block generates a scatter plot of CTD x SG for temperature.

# Choosing some colors 
hmp = 25
colormap_name = 'Reds'
num_colors = 25
cmap = plt.cm.get_cmap(colormap_name)
color_indices = np.linspace(0.3, 1, num_colors)
dc = [cmap(i) for i in color_indices]

legend_times = []

# Select the profiles we want to plot
plot_d = np.arange(2.0,14.4,0.5)

for ind in plot_d:
    # 1. Get the minimum time as an xarray DataArray
    first_time_da = gl.where(gl.dives == ind).dropna("sg_data_point").ctd_time_dt64.mean()
    
    # 2. Extract the underlying NumPy scalar value
    time_scalar = first_time_da.values.item() 
    
    # 3. Check if the value is a valid time (i.e., not Not a Time/NaN)
    # np.isnat checks if a numpy datetime object is NaT
    if pd.notna(time_scalar):
        # Convert the numpy.datetime64/pandas NaT back to a pandas Timestamp for reliable formatting
        first_time_timestamp = pd.to_datetime(time_scalar)
        
        # Format the time (e.g., 'Sep 17, 10:00')
        legend_times.append("SG, "+first_time_timestamp.strftime('%b %d, %H:%M'))
    else:
        # Handle cases where the dive has no associated time (e.g., all times were null)
        legend_times.append('No Data Time')

# Figure begins here! ----------------------------------------
fig, ax = plt.subplots(1, 1, figsize=(10, 10), constrained_layout=True, sharey=True)
    
# Glider Profiles
glider_handles = []
for i, ind in enumerate(range(0,hmp)):
    # Use a label only if it's the main plot for legend handles
    label_text = legend_times[i]
    
    ax.scatter(
        ctd.T,
        ds_bin.temperature_binned[:,ind],
        c=dc[i], 
        zorder=25-ind,
        edgecolor='k',
        lw=2, 
        label=label_text
    )
    glider_line = ax.scatter(
        ctd.T,
        ds_bin.temperature_binned[:,ind],
        c=dc[i], 
        zorder=25-ind
    )
    glider_handles.append(glider_line)

ax.plot([4,18],[4,18],'w',ls='-',lw=3,zorder=100)
ax.plot([4,18],[4,18],'k',ls='--',lw=2,zorder=101)

ax.plot([5,19],[4,18],'.5',ls=':',lw=1, zorder=101)
ax.plot([6,20],[4,18],'.5',ls='-',lw=1, zorder=101)
ax.plot([7,21],[4,18],'.5',ls=':',lw=1, zorder=101)
ax.plot([8,22],[4,18],'.5',ls='-',lw=1, zorder=101)

ax.plot([4,18],[5,19],'.5',ls=':',lw=1, zorder=101)
ax.plot([4,18],[6,20],'.5',ls='-',lw=1, zorder=101)
ax.plot([4,18],[7,21],'.5',ls=':',lw=1, zorder=101)
ax.plot([4,18],[8,22],'.5',ls='-',lw=1, zorder=101)


ax.set_xlim(4,18)
ax.set_ylim(4,18)

# 3. Adding legend to ax (Temperature)
ax.legend(
    glider_handles, 
    legend_times, 
    ncol=2,
    #title='Profile Start Time',
    loc='lower right', 
    bbox_to_anchor=(0.98, 0.02), # Near the bottom-right corner of the plot area
    fontsize=12,
    frameon=False,
    markerfirst=False
    
)

ax.set_xlabel(r"Temperature CTD ($^{\circ}$C)")
ax.set_ylabel(r"Temperature SG ($^{\circ}$C)")

fig.savefig(path_fig+'CTD_cal_Tscatter_Gough_Melktert_Sep2025.png',dpi=300)

# =============================================================================
# Salinity
# =============================================================================

# This block generates a scatter plot of CTD x SG for salinity.

# Choosing some colors 
hmp = 25
colormap_name = 'Blues'
num_colors = 25
cmap = plt.cm.get_cmap(colormap_name)
color_indices = np.linspace(0.3, 1, num_colors)
dc = [cmap(i) for i in color_indices]

legend_times = []
plot_d = np.arange(2.0,14.4,0.5)

for ind in plot_d:
    # 1. Get the minimum time as an xarray DataArray
    first_time_da = gl.where(gl.dives == ind).dropna("sg_data_point").ctd_time_dt64.mean()
    
    # 2. Extract the underlying NumPy scalar value
    time_scalar = first_time_da.values.item() 
    
    # 3. Check if the value is a valid time (i.e., not Not a Time/NaN)
    # np.isnat checks if a numpy datetime object is NaT
    if pd.notna(time_scalar):
        # Convert the numpy.datetime64/pandas NaT back to a pandas Timestamp for reliable formatting
        first_time_timestamp = pd.to_datetime(time_scalar)
        
        # Format the time (e.g., 'Sep 17, 10:00')
        legend_times.append("SG, "+first_time_timestamp.strftime('%b %d, %H:%M'))
    else:
        # Handle cases where the dive has no associated time (e.g., all times were null)
        legend_times.append('No Data Time')

# Figure begins here! ----------------------------------------
fig, ax = plt.subplots(1, 1, figsize=(10, 10), constrained_layout=True, sharey=True)

# Glider Profiles
glider_handles = []
for i, ind in enumerate(range(0,hmp)):
    # Use a label only if it's the main plot for legend handles
    label_text = legend_times[i]
    
    ax.scatter(
        ctd.S,
        ds_bin.salinity_binned[:,ind],
        c=dc[i], 
        zorder=25-ind,
        edgecolor='k',
        lw=2, 
        label=label_text
    )
    glider_line = ax.scatter(
        ctd.S,
        ds_bin.salinity_binned[:,ind],
        c=dc[i], 
        zorder=25-ind#,
        #label=label_text
    )
    glider_handles.append(glider_line)

ax.plot([34.2,35.6],[34.2,35.6],'w',ls='-',lw=3,zorder=100)
ax.plot([34.2,35.6],[34.2,35.6],'k',ls='--',lw=2,zorder=101)

ax.plot([34.3,35.7],[34.2,35.6],'.5',ls=':',lw=1, zorder=101)
ax.plot([34.4,35.8],[34.2,35.6],'.5',ls='-',lw=1, zorder=101)
ax.plot([34.5,35.9],[34.2,35.6],'.5',ls=':',lw=1, zorder=101)
ax.plot([34.6,36],[34.2,35.6],'.5',ls='-',lw=1, zorder=101)

ax.plot([34.2,35.6],[34.3,35.7],'.5',ls=':',lw=1, zorder=101)
ax.plot([34.2,35.6],[34.4,35.8],'.5',ls='-',lw=1, zorder=101)
ax.plot([34.2,35.6],[34.5,35.9],'.5',ls=':',lw=1, zorder=101)
ax.plot([34.2,35.6],[34.6,36],'.5',ls='-',lw=1, zorder=101)

ax.set_xlim(34.4,35.6)
ax.set_ylim(34.4,35.6)

# 3. Adding legend to ax (Salinity)
ax.legend(
    glider_handles, 
    legend_times, 
    ncol=2,
    #title='Profile Start Time',
    loc='lower right', 
    bbox_to_anchor=(1.02, 0.0), # Near the bottom-right corner of the plot area
    fontsize=12,
    frameon=False,
    markerfirst=False
    
)

ax.set_xlabel(r"Salinity CTD (PSU)")
ax.set_ylabel(r"Salinity SG (PSU)")

fig.savefig(path_fig+'CTD_cal_Sscatter_Gough_Melktert_Sep2025.png',dpi=300)

#%%

# =============================================================================
# Linear Regression to calculate the bias and basic statistics
# =============================================================================

# =================
# Temperature

# ---------------------------------------
# Code to perform the linear regression with all 25 first profiles,
# profiles 13 and 13.5 are the first that we have that reach to 1000m 
x = ds_bin.temperature_binned[:, :26].values.ravel()
y =(np.ones_like(ds_bin.temperature_binned[:, :26].values)*ctd.T.values[:,None]).ravel()

mask = np.isfinite(x) & np.isfinite(y)

slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
r2 = r_value**2
rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
bias = np.mean(x[mask] - y[mask])               
x_mean = np.mean(x[mask])
x_std = np.std(x[mask])
y_mean = np.mean(y[mask])
y_std = np.std(y[mask])

print("Using the first 25 TEMPERATURE profiles:\nThe bias is", bias, "\nThe RMSE is", rmse, "\nThe correlation is",r2)
print("--------------------------------------")

# ---------------------------------------
# Code to perform the linear regression with profile number 13 and 13.5,
# between 400 and 900 m 
z1 = 400
z2 = 901

x = ds_bin.temperature_binned[z1:z2, 24:26].values.ravel()
y =(np.ones_like(ds_bin.temperature_binned[z1:z2, 24:26].values)*ctd.T.values[z1:z2,None]).ravel()

mask = np.isfinite(x) & np.isfinite(y)

slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
r2 = r_value**2
rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
bias = np.mean(x[mask] - y[mask])               
x_mean = np.mean(x[mask])
x_std = np.std(x[mask])
y_mean = np.mean(y[mask])
y_std = np.std(y[mask])

print("Using TEMPERATURE profiles 13 and 13.5 between 400 and 900m:\nThe bias is", bias, "\nThe RMSE is", rmse, "\nThe correlation is",r2)
print("--------------------------------------")

# =================
# Salinity

# ---------------------------------------
# Code to perform the linear regression with all 25 first profiles,
# profiles 13 and 13.5 are the first that we have that reach to 1000m 
x = ds_bin.salinity_binned[:, :26].values.ravel()
y =(np.ones_like(ds_bin.salinity_binned[:, :26].values)*ctd.S.values[:,None]).ravel()

mask = np.isfinite(x) & np.isfinite(y)

slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
r2 = r_value**2
rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
bias = np.mean(x[mask] - y[mask])               
x_mean = np.mean(x[mask])
x_std = np.std(x[mask])
y_mean = np.mean(y[mask])
y_std = np.std(y[mask])

print("Using the first 25 SALINITY profiles:\nThe bias for is", bias, "\nThe RMSE is", rmse, "\nThe correlation is",r2)
print("--------------------------------------")

# ---------------------------------------
# Code to perform the linear regression with profile number 13 and 13.5,
# between 400 and 900 m 
z1 = 400
z2 = 901

x = ds_bin.salinity_binned[z1:z2, 24:26].values.ravel()
y =(np.ones_like(ds_bin.salinity_binned[z1:z2, 24:26].values)*ctd.S.values[z1:z2,None]).ravel()

mask = np.isfinite(x) & np.isfinite(y)

slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
r2 = r_value**2
rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
bias = np.mean(x[mask] - y[mask])               
x_mean = np.mean(x[mask])
x_std = np.std(x[mask])
y_mean = np.mean(y[mask])
y_std = np.std(y[mask])

print("Using SALINITY profiles 13 and 13.5 between 400 and 900m:\nThe bias is", bias, "\nThe RMSE is", rmse, "\nThe correlation is",r2)
print("--------------------------------------")

#%%

# =============================================================================
# Chlorophyll
# =============================================================================

colormap_name = 'Greens'
num_colors = 6
cmap = plt.cm.get_cmap(colormap_name)
color_indices = np.linspace(0.3, 1, num_colors)
dc = [cmap(i) for i in color_indices]

# plot_dives is already defined above -- same for every plot
legend_times = []
for ind in plot_dives:
    # 1. Get the minimum time as an xarray DataArray
    first_time_da = gl.where(gl.dives == ind).dropna("sg_data_point").ctd_time_dt64.mean()
    
    # 2. Extract the underlying NumPy scalar value
    time_scalar = first_time_da.values.item() 
    
    # 3. Check if the value is a valid time (i.e., not Not a Time/NaN)
    # np.isnat checks if a numpy datetime object is NaT
    if pd.notna(time_scalar):
        # Convert the numpy.datetime64/pandas NaT back to a pandas Timestamp for reliable formatting
        first_time_timestamp = pd.to_datetime(time_scalar)
        
        # Format the time (e.g., 'Sep 17, 10:00')
        legend_times.append("SG, "+first_time_timestamp.strftime('%b %d, %H:%M'))
    else:
        # Handle cases where the dive has no associated time (e.g., all times were null)
        legend_times.append('No Data Time')

# Combine ship label and time labels for the final legend
all_labels = ['Ship CTD, Sep 17 11:16'] + legend_times
# ----------------------------------------

# Figure begins here! --------------------------------
fig, ax = plt.subplots(1, 1, figsize=(5, 10), constrained_layout=True, sharey=True)

# Function to plot all lines and collect handles
def plot_profiles(current_ax, x_data, y_data, gl_x_col, is_main_plot=False):
    # Ship CTD (Black line) - Only labeled on the main plot for legend collection
    current_ax.plot(x_data, y_data, 'k', lw=4)
    ship_line, = current_ax.plot(x_data, y_data, '0.5', lw=3)
    
    # Glider Profiles
    glider_handles = []
    for i, ind in enumerate(plot_dives):
        # Use a label only if it's the main plot for legend handles
        label_text = legend_times[i] if is_main_plot else None 
        
        glider_line, = current_ax.plot(
            gl[gl_x_col].where(gl.dives == ind), 
            -gl.ctd_pressure.where(gl.dives == ind), 
            c=dc[i], 
            lw=2, 
            label=label_text
        )
        if is_main_plot:
            glider_handles.append(glider_line)
    
    if is_main_plot:
        # Return ship handle and all glider handles
        return [ship_line] + glider_handles
    return None # Return None for inset/salinity plots


# --- ax: Chlorophyll Profile (Main Plot) ---
# Pass is_main_plot=True to collect handles
all_handles = plot_profiles(ax, ctd.flu, -ctd.press, 'wlbb2fl_sig695nm_adjusted', is_main_plot=True)

ax.set_ylim(-1000,0) 
ax.set_xlim(-0.05,1.6)
ax.set_yticks(np.arange(-1000,100,100),labels=np.arange(1000,-100,-100))
#ax.invert_yaxis()
ax.set_xlabel(r'Chlorophyll (mg$\,$m$^{-3}$)')
ax.set_ylabel('Pressure (dbar)')
ax.grid(lw=0.2)

# Adding legend to ax
ax.legend(
    all_handles, 
    all_labels, 
    #title='Profile Start Time',
    loc='lower right', 
    bbox_to_anchor=(0.98, 0.02), # Near the bottom-right corner of the plot area
    fontsize=14,
    frameon=False,
    markerfirst=False
    
)

# --- ax inset: Chl Zoom (0-100m) ---
ax_zoom = ax.inset_axes([0.5, 0.5, 0.5, 0.3]) 
plot_profiles(ax_zoom, ctd.flu, -ctd.press, 'wlbb2fl_sig695nm_adjusted')
ax_zoom.set_ylim(-100, 0)
ax_zoom.set_xlim(-0.05,1.6)
#ax_zoom.invert_yaxis()
ax_zoom.tick_params(labelsize=8)
ax_zoom.grid(lw=0.2)
ax_zoom.set_yticks(np.arange(-100,10,10),labels=np.arange(100,-10,-10))
#ax_zoom.set_yticks([])

ax.indicate_inset_zoom(ax_zoom,edgecolor="k",lw=2)

fig.savefig(path_fig+'/CTD_cal_chl_profiles_Gough_Melktert_Sep2025.png',dpi=300)

#%%

# This block generates a scatter plot of CTD x SG for chlorophyll.

hmp = num_colors = 25
colormap_name = 'Greens'
cmap = plt.cm.get_cmap(colormap_name)
color_indices = np.linspace(0.2, 1, num_colors)
dc = [cmap(i) for i in color_indices]

legend_times = []
plot_d = np.arange(2.0,14.4,0.5)

for ind in plot_d:
    # 1. Get the minimum time as an xarray DataArray
    first_time_da = gl.where(gl.dives == ind).dropna("sg_data_point").ctd_time_dt64.mean()
    
    # 2. Extract the underlying NumPy scalar value
    time_scalar = first_time_da.values.item() 
    
    # 3. Check if the value is a valid time (i.e., not Not a Time/NaN)
    # np.isnat checks if a numpy datetime object is NaT
    if pd.notna(time_scalar):
        # Convert the numpy.datetime64/pandas NaT back to a pandas Timestamp for reliable formatting
        first_time_timestamp = pd.to_datetime(time_scalar)
        
        # Format the time (e.g., 'Sep 17, 10:00')
        legend_times.append("SG, "+first_time_timestamp.strftime('%b %d, %H:%M'))
    else:
        # Handle cases where the dive has no associated time (e.g., all times were null)
        legend_times.append('No Data Time')

fig, ax = plt.subplots(1, 1, figsize=(10, 10), constrained_layout=True, sharey=True)
    
# Glider Profiles 
glider_handles = []
for i, ind in enumerate(range(0,hmp)):
    # Use a label only if it's the main plot for legend handles
    label_text = legend_times[i]
    
    nans = (np.isnan(ctd.flu.values) | np.isnan(ds_bin.wlbb2fl_sig695nm_adjusted_binned[:,ind]))
    smoo_chl_ctd = smoo(ctd.flu.values[~nans],10)
    smoo_chl_gl = smoo(ds_bin.wlbb2fl_sig695nm_adjusted_binned[:,ind][~nans],10)
    
    ax.scatter(
        np.log10(smoo_chl_ctd),
        np.log10(smoo_chl_gl),
        c=dc[i], 
        zorder=25-ind,
        edgecolor='k',
        lw=2, 
        label=label_text
    )
    glider_line = ax.scatter(
        np.log10(smoo_chl_ctd),
        np.log10(smoo_chl_gl),
        c=dc[i], 
        zorder=25-ind
    )
    glider_handles.append(glider_line)

ax.plot([-2,.5],[-2,.5],'w',ls='-',lw=3,zorder=100)
ax.plot([-2,.5],[-2,.5],'k',ls='--',lw=2,zorder=101)

ax.plot([-2.5,0],[-2,.5],'.5',ls=':',lw=1, zorder=101)
ax.plot([-3,-0.5],[-2,.5],'.5',ls='-',lw=1, zorder=101)
ax.plot([-3.5,-1],[-2,.5],'.5',ls=':',lw=1, zorder=101)
ax.plot([-4,-1.5],[-2,.5],'.5',ls='-',lw=1, zorder=101)

ax.plot([-2,.5],[-2.5,0],'.5',ls=':',lw=1, zorder=101)
ax.plot([-2,.5],[-3,-0.5],'.5',ls='-',lw=1, zorder=101)
ax.plot([-2,.5],[-3.5,-1],'.5',ls=':',lw=1, zorder=101)
ax.plot([-2,.5],[-4,-1.5],'.5',ls='-',lw=1, zorder=101)

ax.set_xlim(-2,0.5)
ax.set_ylim(-2,0.5)

# Adding legend to ax (Chlorophyll)
ax.legend(
    glider_handles, 
    legend_times, 
    ncol=1,
    loc='lower right', 
    bbox_to_anchor=(0.26, 0.3), # Near the bottom-right corner of the plot area
    fontsize=12,
    frameon=False,
    markerfirst=False
    
)

ax.set_xlabel(r"Chlorophyll CTD (mg$\,$m$^{-3}$)")
ax.set_ylabel(r"Chlorophyll SG (mg$\,$m$^{-3}$)")

ax.set_yticks([0,-1,-2])
ax.set_xticks([0,-1,-2])

ax.set_yticklabels([r"10$^{-2}$",r"10$^{-1}$",r"10$^{0}$"])
ax.set_xticklabels([r"10$^{-2}$",r"10$^{-1}$",r"10$^{0}$"])

fig.savefig(path_fig+'CTD_cal_chl_scatter_Gough_Melktert_Sep2025.png',dpi=300)

#%%

# =============================================================================
# Linear Regression to calculate the bias and basic statistics
# =============================================================================

# =================
# Chlorophyll

# ---------------------------------------
# Code to perform the linear regression with all 25 first profiles,
# profiles 13 and 13.5 are the first that we have that reach to 1000m 
x = ds_bin.wlbb2fl_sig695nm_adjusted_binned[:, :26].values.ravel()
y =(np.ones_like(ds_bin.wlbb2fl_sig695nm_adjusted_binned[:, :26].values)*ctd.flu.values[:,None]).ravel()

mask = np.isfinite(x) & np.isfinite(y)

slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
r2 = r_value**2
rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
bias = np.mean(x[mask] - y[mask])               
x_mean = np.mean(x[mask])
x_std = np.std(x[mask])
y_mean = np.mean(y[mask])
y_std = np.std(y[mask])

print("Using the first 25 CHLOROPHYL profiles:\nThe bias is", bias, "\nThe RMSE is", rmse, "\nThe correlation is",r2)
print("--------------------------------------")

# ---------------------------------------
# Code to perform the linear regression with profile number 13 and 13.5,
# between 400 and 900 m 
z1 = 400
z2 = 901

x = ds_bin.wlbb2fl_sig695nm_adjusted_binned[z1:z2, :26].values.ravel()
y =(np.ones_like(ds_bin.wlbb2fl_sig695nm_adjusted_binned[z1:z2, :26].values)*ctd.flu.values[z1:z2,None]).ravel()

mask = np.isfinite(x) & np.isfinite(y)

slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
r2 = r_value**2
rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
bias = np.mean(x[mask] - y[mask])               
x_mean = np.mean(x[mask])
x_std = np.std(x[mask])
y_mean = np.mean(y[mask])
y_std = np.std(y[mask])

print("Using CHLOROPHYL profiles 13 and 13.5 between 400 and 900m:\nThe bias is", bias, "\nThe RMSE is", rmse, "\nThe correlation is",r2)
print("--------------------------------------")

#%%
# =============================================================================
# Stats -- all profiles
# =============================================================================

# Example:
# dic = [
#     ("Dive 2", ds_bin.sel(profile=2), ctd,["temperature_binned","salinity_binned"]),
#     ("Dive 3", ds_bin.sel(profile=3), ctd,["temperature_binned","salinity_binned"]),
#     ("Dive 8", ds_bin.sel(profile=8), ctd,["temperature_binned","salinity_binned"]),
#     ("Dive 13", ds_bin.sel(profile=13), ctd,["temperature_binned","salinity_binned"]),
#     ("Dive 13.5", ds_bin.sel(profile=13.5), ctd,["temperature_binned","salinity_binned"]),
# ]

dic = []

# Calculate the bias for all profiles, no depth constraint.
# If you want to constrain by depth, replace ds_bin.sel(profile=i) by ds_bin.sel(profile=i,pressure=slice(400,900))
# and ctd by ctd.sel(index=slice(400,900))

for i in ds_bin.profile.values:
    dic.append((f"Dive {i}", ds_bin.sel(profile=i), ctd,["temperature_binned","salinity_binned"]))

ship_var = ["T","S"]

stat_dfs = {}

all_bias_T = []
all_bias_S = []

for DIVE, dsx, dsy, comp_vars in dic:
    stats = []
    for x_var, y_var in zip(comp_vars, ship_var):
        if x_var in dsx and y_var in dsy:
            x = dsx[x_var].values
            y = dsy[y_var].values
            mask = np.isfinite(x) & np.isfinite(y)
            n = mask.sum()
            if n > 0:
                slope, intercept, r_value, p_value, std_err = linregress(x[mask], y[mask])
                r2 = r_value**2
                rmse = np.sqrt(np.mean((x[mask] - y[mask])**2))
                bias = np.mean(x[mask] - y[mask])               
                x_mean = np.mean(x[mask])
                x_std = np.std(x[mask])
                y_mean = np.mean(y[mask])
                y_std = np.std(y[mask])
            else:
                slope = intercept = r2 = p_value = rmse = bias = np.nan
                x_mean = x_std = y_mean = y_std = np.nan
            stats.append({
                "variable": x_var,
                "n": n,
                "r2": r2,
                "slope": slope,
                "intercept": intercept,
                "p": p_value,
                "rmse": rmse,
                "bias": bias,})
            
            if x_var =='salinity_binned':
                all_bias_S.append(bias)
            else: all_bias_T.append(bias)
        
    stat_dfs[DIVE] = pd.DataFrame(stats)

#%%

# =============================================================================
# Calculating Earth-relative distances based on lat and lon and plotting the biases
# =============================================================================

lat = np.nanmean(ds_bin.longitude_binned.values,axis=0)
lon = np.nanmean(ds_bin.latitude_binned.values,axis=0)

nans = (np.isnan(lat) & np.isnan(lon))

# Distance vector
dist = np.cumsum(gsw.distance(lon[~nans],lat[~nans]))

# Time vector
t=[]
for i in np.arange(1,261,0.5): # see what's going on in here...
    t.append(np.nanmean(gl.where(gl.dives==i).ctd_time.values))

t = np.array(t)
nans = np.isnan(t)
time = get_date(t[~nans])

#-----------------

ind = 26 # profiles that should be plotted in colour
figure_name = 'biases_TS_CTD_cal_Gough_Melktert_Sep2025.png'
# Figure begins here! ----------------------
fig, [ax, bx] = plt.subplots(2, 1, figsize=(12, 10), constrained_layout=True, sharex=True)

ax.plot(dist/1e3,all_bias_T[3:],'k',lw=1)
ax.plot(dist[:ind]/1e3,all_bias_T[3:ind+3],'r',lw=3,marker='.',markersize=10)
ax.set_ylabel("Temperature bias ($^{\circ}$C)")
ax.set_ylim([-2.5,0.5])
ax.set_xticks(np.arange(0,1050,100))
ax.grid(lw=0.2)

bx.plot(dist/1e3,all_bias_S[3:],'k',lw=1)
bx.plot(dist[:ind]/1e3,all_bias_S[3:ind+3],'b',lw=3,marker='.',markersize=10)
bx.set_ylabel("Salinity bias (PSU)")
bx.set_ylim([-0.25,0.05])
bx.set_xlim([0,1040])
bx.set_xticks(np.arange(0,1050,100))
bx.set_xticklabels(["0", "100",  "200",  "300",  "400",  "500",  "600",  "700",  "800",  "900", "1000"])
bx.set_xlabel("Distance (km)")
bx.grid(lw=0.2)

fig.savefig(path_fig+figure_name,dpi=300)