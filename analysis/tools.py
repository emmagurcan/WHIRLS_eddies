import matplotlib.pyplot as plt
import h5py
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os as os 
import numpy as np
import cmocean
import netCDF4 as nc
from cartopy.feature import LAND
from cartopy.feature import ShapelyFeature
import xarray as xr

import geopandas as gpd
import os as os 
import scipy
import cfgrib
# import fonctions as funcs

#For TOEddies
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
from pyproj import Geod
import numpy as np
import matplotlib.lines as mlines
from matplotlib.patches import FancyArrowPatch
def dm_to_dd(degrees, minutes, flag=False):
    dd = degrees + (minutes / 60)
    if flag:
        dd = -dd
    return dd

def haversine_nm(lon1, lat1, lon2, lat2):
    R = 6371.0  # km
    deg2rad = np.pi / 180.0

    dlon = (lon2 - lon1) * deg2rad
    dlat = (lat2 - lat1) * deg2rad

    a = (np.sin(dlat/2)**2 +
         np.cos(lat1 * deg2rad) * np.cos(lat2 * deg2rad) *
         np.sin(dlon/2)**2)

    c = 2 * np.arcsin(np.sqrt(a))
    km = R * c
    nm = km * 0.539957  # convert km to nautical miles
    return nm
def format_time(h):
    hours = int(h)
    minutes = int(round((h - hours) * 60))
    return hours, minutes

def get_time_mask(month, day,new_adcp, t_adcp, dm = True):
    mask = (t_adcp.dt.month == month) & (t_adcp.dt.day == day)
    new_adcp = new_adcp.isel(time=mask)
    if dm:
        depth_mask = (new_adcp.depth < -50) & (new_adcp.depth > -100)
        new_adcp = new_adcp.isel(depth=depth_mask)
    return new_adcp

def get_depth_mask(new_adcp):
    depth_mask = (new_adcp.depth < -50) & (new_adcp.depth > -100)
    new_adcp = new_adcp.isel(depth=depth_mask)
    return new_adcp

geod = Geod(ellps="WGS84")

def along_track_distance(lon, lat):
    """
    Returns cumulative distance in meters along a track.
    """
    lon = np.asarray(lon)
    lat = np.asarray(lat)

    # distance between consecutive points
    _, _, dist = geod.inv(lon[:-1], lat[:-1], lon[1:], lat[1:])

    # cumulative sum (prepend 0 for starting point)
    return np.concatenate([[0], np.cumsum(dist)])


def load_eddies(dir_eddies):
    """Load cyclonic and anticyclonic eddy shapefiles."""

    anti = {}
    cyc = {}

    for fname in dir_eddies:

        if "anti" in str(fname):
            target = anti
        else:
            target = cyc

        if "maxspeed" in str(fname):
            target["maxspeed"] = gpd.read_file(fname)

        elif "outcontour" in str(fname):
            target["outcontour"] = gpd.read_file(fname)

        elif "centers" in str(fname):
            target["centers"] = gpd.read_file(fname)

    return anti, cyc

def load_satellite(adt_path):

    # ds = xr.open_dataset(adt_daily_dir + file_adt)
    ds = xr.open_dataset(adt_path)
    if "time" in ds:
        ds = ds.drop_vars("time")
    ds.close()
    return {
        "dataset": ds,
        "lon": ds.longitude,
        "lat": ds.latitude,
        "adt": ds.adt,
        "ugos": np.squeeze(ds.ugos),
        "vgos": np.squeeze(ds.vgos)
    }
from pathlib import Path
from datetime import datetime

def get_closest_sst_file(date, sst_dir):
    """
    Return the SST file corresponding to `date`, or the closest available date.

    Parameters
    ----------
    date : str
        Date in the format "YYYY-MM-DD".
    sst_dir : str or Path
        Directory containing the SST files.

    Returns
    -------
    str
        Full path to the selected SST file.
    """

    sst_dir = Path(sst_dir)

    target_date = datetime.strptime(date, "%Y-%m-%d")

    # Exact filename
    exact_file = sst_dir / f"BT_masked_Meteosat_{target_date.strftime('%Y%m%d')}.nc"
    print(exact_file)
    if exact_file.exists():
        return str(exact_file)

    # Search all SST files
    files = sorted(sst_dir.glob("BT_masked_Meteosat_*.nc"))

    if len(files) == 0:
        raise FileNotFoundError(f"No SST files found in {sst_dir}")

    available_dates = []
    for f in files:
        try:
            date_str = f.stem.split("_")[-1]   # e.g. 20260710
            file_date = datetime.strptime(date_str, "%Y%m%d")
            available_dates.append((file_date, f))
        except ValueError:
            continue

    if len(available_dates) == 0:
        raise FileNotFoundError("No valid SST files found.")

    # Find the closest date
    closest_date, closest_file = min(
        available_dates,
        key=lambda x: abs((x[0] - target_date).days)
    )

    print(f"No SST file for {date}. Using {closest_date.strftime('%Y-%m-%d')} instead.")

    return str(closest_file)

def load_satellite(adt_path):

    # ds = xr.open_dataset(adt_daily_dir + file_adt)
    ds = xr.open_dataset(adt_path)
    if "time" in ds:
        ds = ds.drop_vars("time")
    ds.close()
    return {
        "dataset": ds,
        "lon": ds.longitude,
        "lat": ds.latitude,
        "adt": ds.adt,
        "ugos": np.squeeze(ds.ugos),
        "vgos": np.squeeze(ds.vgos)
    }


def average_adcp_velocity(curr_day):

    u = curr_day.uvel_adcp.values
    v = curr_day.vvel_adcp.values

    mask = (
        ~np.isnan(u).any(axis=1)
        &
        ~np.isnan(v).any(axis=1)
    )

    u_mean = np.full(u.shape[0], np.nan)
    v_mean = np.full(v.shape[0], np.nan)

    u_mean[mask] = u[mask].mean(axis=1)
    v_mean[mask] = v[mask].mean(axis=1)

    return u_mean, v_mean

def cut_adcp(adcp, x1_bound, x2_bound, x_dist):
    x_km = x_dist / 1000

    mask = (x_km <= x1_bound) | (x_km >= x2_bound)

    return adcp.isel(time=mask), x_dist[mask]

def shift_results(results, x_start, offset):

    shifted = []

    if results is not None:
        for r in results:
            if type(r) is not str:
                r_new = r.copy()
            else:
                r_new = r

            # move from original x coordinate
            # to glued x coordinate
            r_new["x_dist"] = (
                r["x_dist"] - x_start + offset
            )

            shifted.append(r_new)

    return shifted

def velocity_magnitude(u, v):
    return np.sqrt(u**2 + v**2)

def interpolate_satellite(curr_day, ugos, vgos):

    u = ugos.interp(
        longitude=curr_day.lon_adcp,
        latitude=curr_day.lat_adcp,
        method="nearest"
    )

    v = vgos.interp(
        longitude=curr_day.lon_adcp,
        latitude=curr_day.lat_adcp,
        method="nearest"
    )

    return u, v

def plot_eddies(ax, eddies, color, crs):

    ms = eddies["maxspeed"]
    oc = eddies["outcontour"]
    cent = eddies["centers"]

    for i in range(len(cent)):

        if ms.geometry[i] is None:
            continue

        x, y = ms.geometry[i].exterior.xy
        ax.plot(x, y, color=color, transform=crs, zorder=20)

        x, y = oc.geometry[i].exterior.xy
        ax.plot(
            x,
            y,
            "--",
            color=color,
            transform=crs,
            zorder=20
        )

        ax.scatter(
            cent.geometry[i].x,
            cent.geometry[i].y,
            s=10,
            color=color,
            transform=crs,
            zorder=20,
        )
        # ax.text(
        #     cent.geometry[i].x + 0.03,   # longitude offset
        #     cent.geometry[i].y + 0.03,   # latitude offset
        #     str(cent.eddy_index[i]),
        #     transform=crs,
        #     fontsize=10,
        #     color=color,
        #     zorder=200,
        #     bbox=dict(
        #         boxstyle="round,pad=0.2",
        #         facecolor="white",
        #         edgecolor=color,
        #         alpha=0.8,
        #     ),
        # )


def plot_along_track(
    x_dist,
    adcp_mag,
    sat_mag,
    month,
    day
):

    plt.figure()

    plt.plot(
        -x_dist/1000,
        adcp_mag,
        color="blueviolet",
        label="ADCP"
    )

    plt.plot(
        -x_dist/1000,
        sat_mag,
        color="orange",
        label="Satellite"
    )

    plt.legend()
    plt.xlabel("Along-track distance (km)")
    plt.ylabel("Velocity (m/s)")
    plt.title(f"{month}-{day}")

def plot_adcp_section(curr_day, adcp_mag, x_dist, title, vmin, vmax, cmap="viridis", cyc = True, results_oc=None, results_ms=None, eddy_center_toeddies=None, eddy_center_nencioli=None, savefig=None):

    fig, ax = plt.subplots(figsize=(9,4))

    c = ax.pcolormesh(
        x_dist,
        curr_day.depth.values,
        adcp_mag.T,
        vmin = vmin,
        vmax = vmax,
        cmap=cmap
    )
    if cyc:
        col = 'blue'
    else:
        col = 'red'
    if results_oc is not None:
        for result in results_oc:
            ax.axvline(result["x_dist"], color=col, linestyle='--')
    if results_ms is not None:
        for result in results_ms:
            ax.axvline(result["x_dist"], color=col, linestyle='solid')
    if eddy_center_toeddies is not None:
        ax.axvline(eddy_center_toeddies["x_dist"], color=col, linestyle='dashdot')

    if eddy_center_nencioli is not None:
        ax.axvline(eddy_center_nencioli, color=col, linestyle='dashdot')
    ax.invert_xaxis()
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label("Velocity [m/s]")
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color=col, linestyle='--', label="outer contour"),
        Line2D([0], [0], color=col, linestyle='solid', label='max speed'),
        Line2D([0], [0], color=col, linestyle='dashdot', label='center'),
    ]

    leg = ax.legend(
        handles=legend_handles,
        loc="lower left",
        fontsize="small",
        frameon=True,
    )
    ax.set_xlabel("Distance [km]")
    ax.set_ylabel("Depth [m]")
    ax.set_title(title)
    if savefig is not None:
        plt.savefig(savefig)

    return fig, ax
def setup_map(ax, extent, title, crs, ft):
    """Apply common formatting to Cartopy maps."""

    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_title(title)

    ax.add_feature(
        cfeature.LAND,
        linewidth=0.5,
        edgecolor="black",
        facecolor="gray",
        zorder=4
    )

    ax.coastlines(
        resolution="10m",
        color="k",
        zorder=5
    )

    grid = ax.gridlines(
        draw_labels=True,
        linewidth=1,
        color="black",
        alpha=0.2,
        linestyle="--",
        zorder=6
    )

    grid.top_labels = False
    grid.right_labels = False
    grid.left_labels = True
    grid.xlabel_style = {"size": ft}
    grid.ylabel_style = {"size": ft}

def plot_adcp_satellite_map(
    curr_day,
    sat,
    sat_mag,
    u_adcp,
    v_adcp,
    ugos_track,
    vgos_track,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    longs_bato,
    lat_bato,
    x_dist,
):
    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    c = ax.contourf(
        sat["lon"],
        sat["lat"],
        sat_mag,
        levels=np.arange(0, 2.5, 0.1),
        cmap="Spectral_r",
        transform=crs,
    )

    ax.scatter(
        longs_bato[::50],
        lat_bato[::50],
        s=1,
        color="black",
        transform=crs,
        zorder=15,
    )

    step = 50

    # =========================
    # ADCP vectors
    # =========================
    q_adcp = ax.quiver(
        curr_day.lon_adcp.values[::step],
        curr_day.lat_adcp.values[::step],
        u_adcp[::step],
        v_adcp[::step],
        transform=crs,
        scale=3,
        width=0.005,
        color="blueviolet",
        zorder=100,
    )

    # =========================
    # Satellite vectors
    # =========================
    q_sat = ax.quiver(
        curr_day.lon_adcp.values[::step],
        curr_day.lat_adcp.values[::step],
        ugos_track.values[::step],
        vgos_track.values[::step],
        transform=crs,
        scale=3,
        width=0.005,
        color="orange",
        zorder=100,
    )

    # =========================
    # 🔥 DISTANCE LABELS (NEW)
    # =========================
    target_dists = np.arange(50, np.nanmax(x_dist)/1000, 50)  # km

    x_km = x_dist / 1000
    lon = curr_day.lon_adcp.values
    lat = curr_day.lat_adcp.values

    for d in target_dists:
        idx = np.argmin(np.abs(x_km - d))

        ax.text(
            lon[idx],
            lat[idx],
            f"{int(d)} km",
            transform=crs,
            fontsize=3,
            color="black",
            ha="left",
            va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="black",
                linewidth=0.2,
                alpha=0.8,
            ),
            zorder=120,
        )

    # =========================
    # Eddies
    # =========================
    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    # =========================
    # Colorbar
    # =========================
    cb = fig.colorbar(c, ax=ax, shrink=0.6)
    cb.set_label("Geostrophic Velocity [m/s]")

    # =========================
    # Reference arrow
    # =========================
    ax.quiverkey(
        q_adcp,
        X=1.1,
        Y=0.9,
        U=0.2,
        label="0.2 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color="black"
    )


    # =========================
    # Legend
    # =========================
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color="blueviolet", lw=2, label="ADCP"),
        Line2D([0], [0], color="orange", lw=2, label="Satellite"),
    ]

    leg = ax.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize="small",
        frameon=True,
    )
    leg.set_zorder(110)

    return fig, ax
# def plot_adcp_satellite_map(
#     curr_day,
#     sat,
#     sat_mag,
#     u_adcp,
#     v_adcp,
#     ugos_track,
#     vgos_track,
#     anti,
#     cyc,
#     extent,
#     title,
#     proj,
#     crs,
#     ft,
#     longs_bato,
#     lat_bato,
# ):
#     fig = plt.figure(dpi=200)
#     ax = fig.add_subplot(111, projection=proj)

#     setup_map(ax, extent, title, crs, ft)

#     c = ax.contourf(
#         sat["lon"],
#         sat["lat"],
#         sat_mag,
#         levels=np.arange(0, 2.5, 0.1),
#         cmap="Spectral_r",
#         transform=crs,
#     )

#     ax.scatter(
#         longs_bato[::50],
#         lat_bato[::50],
#         s=1,
#         color="black",
#         transform=crs,
#         zorder=15,
#     )

#     step = 50

#     # =========================
#     # ADCP vectors
#     # =========================
#     q_adcp = ax.quiver(
#         curr_day.lon_adcp.values[::step],
#         curr_day.lat_adcp.values[::step],
#         u_adcp[::step],
#         v_adcp[::step],
#         transform=crs,
#         scale=3,
#         width=0.005,
#         color="blueviolet",
#         zorder=100,
#     )

#     # =========================
#     # Satellite vectors
#     # =========================
#     q_sat = ax.quiver(
#         curr_day.lon_adcp.values[::step],
#         curr_day.lat_adcp.values[::step],
#         ugos_track.values[::step],
#         vgos_track.values[::step],
#         transform=crs,
#         scale=3,
#         width=0.005,
#         color="orange",
#         zorder=100,
#     )

#     # =========================
#     # Eddies
#     # =========================
#     plot_eddies(ax, anti, "red", crs)
#     plot_eddies(ax, cyc, "blue", crs)

#     # =========================
#     # Colorbar
#     # =========================
#     cb = fig.colorbar(c, ax=ax, shrink=0.6)
#     cb.set_label("Geostrophic Velocity [m/s]")

#     # =========================
#     # Reference arrow (0.5 m/s)
#     # =========================
#     # =========================
#     # Reference arrow (0.5 m/s)
#     # =========================
#     ax.quiverkey(
#         q_adcp,   # or q_sat
#         X=1.1,
#         Y=0.9,
#         U=0.2,
#         label="0.2 m/s",
#         fontproperties={"size": 8},
#         labelpos="E",
#         color="black")
        
#     # =========================
#     # LEGEND (manual handles)
#     # =========================
#     from matplotlib.lines import Line2D

#     legend_handles = [
#         Line2D([0], [0], color="blueviolet", lw=2, label="ADCP"),
#         Line2D([0], [0], color="orange", lw=2, label="Satellite"),
#     ]

#     leg = ax.legend(
#         handles=legend_handles,
#         loc="upper right",
#         fontsize="small",
#         frameon=True,
#     )
#     leg.set_zorder(110) 
#     return fig, ax
def plot_adt_map(
    sat,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    vmin,
    vmax,
    lon_subset,
    lat_subset,
    depth_subset,
    longs_bato,
    lat_bato,
    curr_day,
    x_dist,
    dx=50
):

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    bathy = ax.contour(
        lon_subset,
        lat_subset,
        depth_subset,
        colors="grey",
        linewidths=0.8,
        transform=crs,
    )

    ax.clabel(bathy, fontsize=4)

    c = ax.contourf(
        sat["lon"],
        sat["lat"],
        sat["adt"][0],
        levels=np.arange(vmin, vmax, 0.01),
        cmap="RdBu_r",
        transform=crs,
    )

    lon2d, lat2d = np.meshgrid(
        sat["lon"].values,
        sat["lat"].values
    )

    step = 2

    Q = ax.quiver(
        lon2d[::step, ::step],
        lat2d[::step, ::step],
        sat["ugos"].values[::step, ::step],
        sat["vgos"].values[::step, ::step],
        transform=crs,
        scale=18,
    )
    target_dists = np.arange(np.nanmin(x_dist)/1000, np.nanmax(x_dist)/1000, dx)  # km

    x_km = x_dist / 1000
    lon = curr_day.lon_adcp.values
    lat = curr_day.lat_adcp.values

    for d in target_dists:
        idx = np.argmin(np.abs(x_km - d))
        t_curr = np.datetime_as_string(curr_day.time[idx].values,  unit='D')
        ax.text(
            lon[idx],
            lat[idx],
            f"{int(d)} km \n {t_curr}",
            transform=crs,
            fontsize=3,
            color="black",
            ha="left",
            va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="black",
                linewidth=0.2,
                alpha=0.8,
            ),
            zorder=120,
        )

    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    ax.scatter(
        longs_bato[::50],
        lat_bato[::50],
        s=0.5,
        color="green",
        transform=crs,
        zorder=90,
    )

    ax.quiverkey(
        Q,
        X=1.08,
        Y=0.92,
        U=0.5,
        label="0.5 m/s",
        labelpos="E",
    )

    cb = fig.colorbar(c, ax=ax, shrink=0.6)
    cb.set_label("ADT [m]")

    return fig, ax


def plot_alongtrack(x_dist, adcp_mag, sat_mag, title):

    fig, ax = plt.subplots(figsize=(8,4))

    ax.plot(
        -x_dist / 1000,
        adcp_mag,
        color="blueviolet",
        label="ADCP",
    )

    ax.plot(
        -x_dist / 1000,
        sat_mag,
        color="orange",
        label="Satellite",
    )

    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Velocity (m/s)")
    ax.set_title(title)
    ax.legend()

    return fig, ax
def plot_sst_map(
    sat,
    lon,
    lat,
    sst,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    vmin,
    vmax,
    lon_subset,
    lat_subset,
    depth_subset,
    longs_bato,
    lat_bato,
    curr_day,
    x_dist,
    dx=50
):

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    bathy = ax.contour(
        lon_subset,
        lat_subset,
        depth_subset,
        colors="grey",
        linewidths=0.8,
        transform=crs,
    )

    ax.clabel(bathy, fontsize=4)

    c = ax.pcolormesh(
        lon,
        lat,
        sst[:,:], vmin=12, vmax=15, cmap=cmocean.cm.thermal,transform=crs
    )

    lon2d, lat2d = np.meshgrid(
        sat["lon"].values,
        sat["lat"].values
    )

    step = 2

    Q = ax.quiver(
        lon2d[::step, ::step],
        lat2d[::step, ::step],
        sat["ugos"].values[::step, ::step],
        sat["vgos"].values[::step, ::step],
        transform=crs,
        scale=18,
    )


    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    ax.scatter(
        longs_bato[::50],
        lat_bato[::50],
        s=0.5,
        color="green",
        transform=crs,
        zorder=90,
    )

    ax.quiverkey(
        Q,
        X=1.08,
        Y=0.94,
        U=0.5,
        label="0.5 m/s",
        labelpos="E",
    )

    cb = fig.colorbar(c, ax=ax, shrink=0.6, extend='both')
    cb.set_label("ADT [m]")

    return fig, ax
def plot_rel_vort_map(
    curr_day,
    u_adcp,
    v_adcp,
    ugos_track,
    vgos_track,
    sat,
    rel_vort,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    vmin,
    vmax,
    lon_subset,
    lat_subset,
    depth_subset,
    longs_bato,
    lat_bato,
    cbar_lab,
    deg=0.125,
    cmap="RdBu_r",
    extend=True
):

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    bathy = ax.contour(
        lon_subset,
        lat_subset,
        depth_subset,
        colors="grey",
        linewidths=0.8,
        transform=crs,
    )

    ax.clabel(bathy, fontsize=4)
    levels = np.linspace(vmin, vmax, 51)
    if extend:
        e = 'both'
    else:
        e = 'min'
    c = ax.contourf(
        sat["lon"],
        sat["lat"],
        rel_vort,
        levels=levels,
        cmap=cmap,
        transform=crs,
        extend=e
    )

    
    lon = curr_day.lon_adcp.values
    lat = curr_day.lat_adcp.values

    u = u_adcp
    v = v_adcp

    # round to nearest 0.1 degree grid
    lon_bin = np.round(lon / deg) * deg
    lat_bin = np.round(lat / deg) * deg

    # keep only unique grid points
    mask = np.ones(lon.shape, dtype=bool)

    seen = set()

    for i in range(lon.size):
        key = (lon_bin[i], lat_bin[i])
        if key in seen:
            mask[i] = False
        else:
            seen.add(key)
    q_adcp = ax.quiver(
        lon[mask],
        lat[mask],
        u[mask],
        v[mask],
        transform=crs,
        scale=8,
        width=0.005,
        color="blueviolet",
        zorder=100,
    )

    lon2d, lat2d = np.meshgrid(
        sat["lon"].values,
        sat["lat"].values
    )

    step = 2

    Q = ax.quiver(
        lon2d[::step, ::step],
        lat2d[::step, ::step],
        sat["ugos"].values[::step, ::step],
        sat["vgos"].values[::step, ::step],
        transform=crs,
        scale=18,
    )

    # =========================
    # Satellite vectors
    # =========================
    q_sat = ax.quiver(
        curr_day.lon_adcp.values[mask],
        curr_day.lat_adcp.values[mask],
        ugos_track.values[mask],
        vgos_track.values[mask],
        transform=crs,
        scale=8,
        width=0.005,
        color="orange",
        zorder=100,
    )

    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    ax.scatter(
        longs_bato[::50],
        lat_bato[::50],
        s=0.5,
        color="black",
        transform=crs,
        zorder=90,
    )
    ax.quiverkey(
        q_adcp,
        X=1.1,
        Y=0.9,
        U=0.5,
        label="0.5 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color="orange"
    )
    ax.quiverkey(
        q_sat,
        X=1.1,
        Y=0.86,
        U=0.5,
        label="0.5 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color="blueviolet"
    )
    ax.quiverkey(
        Q,
        X=1.1,
        Y=0.94,
        U=0.5,
        label="0.5 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color = "black"
    )


    cb = fig.colorbar(c, ax=ax, shrink=0.6)

    ticks = np.arange(vmin, vmax + 0.1, 0.1)
    if "Total Strain" in title:
        ticks = np.arange(vmin, vmax, 0.1)

    cb.set_ticks(ticks)
    cb.set_label(f"{cbar_lab} [s-1]")
    return fig, ax

    
def plot_adcp_satellite_map_no_xdist(
    curr_day,
    sat,
    sat_mag,
    u_adcp,
    v_adcp,
    ugos_track,
    vgos_track,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    longs_bato,
    lat_bato,
    deg=0.2
):
    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    c = ax.contourf(
        sat["lon"],
        sat["lat"],
        sat_mag,
        levels=np.arange(0, 2.5, 0.05),
        cmap="Spectral_r",
        transform=crs,
    )

    ax.scatter(
        longs_bato[::50],
        lat_bato[::50],
        s=1,
        color="black",
        transform=crs,
        zorder=15,
    )

    # =========================
    # ADCP vectors
    # =========================

    lon = curr_day.lon_adcp.values
    lat = curr_day.lat_adcp.values

    u = u_adcp
    v = v_adcp

    # round to nearest 0.1 degree grid
    lon_bin = np.round(lon / 0.2) * 0.2
    lat_bin = np.round(lat / 0.2) * 0.2

    # keep only unique grid points
    mask = np.ones(lon.shape, dtype=bool)

    seen = set()

    for i in range(lon.size):
        key = (lon_bin[i], lat_bin[i])
        if key in seen:
            mask[i] = False
        else:
            seen.add(key)
    q_adcp = ax.quiver(
        lon[mask],
        lat[mask],
        u[mask],
        v[mask],
        transform=crs,
        scale=5,
        width=0.005,
        color="blueviolet",
        zorder=100,
    )
    step = 2
    lon2d, lat2d = np.meshgrid(
        sat["lon"].values,
        sat["lat"].values
    )
    Q = ax.quiver(
        lon2d[::step, ::step],
        lat2d[::step, ::step],
        sat["ugos"].values[::step, ::step],
        sat["vgos"].values[::step, ::step],
        transform=crs,
        scale=18,
    )

    # =========================
    # Satellite vectors
    # =========================
    q_sat = ax.quiver(
        curr_day.lon_adcp.values[mask],
        curr_day.lat_adcp.values[mask],
        ugos_track.values[mask],
        vgos_track.values[mask],
        transform=crs,
        scale=5,
        width=0.005,
        color="orange",
        zorder=100,
    )

    # =========================
    # Eddies
    # =========================
    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    # =========================
    # Colorbar
    # =========================
    cb = fig.colorbar(c, ax=ax, shrink=0.6)
    cb.set_label("Geostrophic Velocity [m/s]")

    # =========================
    # Reference arrow
    # =========================
    ax.quiverkey(
        q_adcp,
        X=1.1,
        Y=0.9,
        U=0.5,
        label="0.5 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color="orange"
    )
    ax.quiverkey(
        q_sat,
        X=1.1,
        Y=0.86,
        U=0.5,
        label="0.5 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color="blueviolet"
    )
    
    ax.quiverkey(
        Q,
        X=1.1,
        Y=0.94,
        U=0.5,
        label="0.5 m/s",
        labelpos="E",
        fontproperties={"size": 8},
        color = "black"
    )


    # =========================
    # Legend
    # =========================
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color="blueviolet", lw=2, label="ADCP"),
        Line2D([0], [0], color="orange", lw=2, label="Satellite"),
    ]

    leg = ax.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize="small",
        frameon=True,
    )
    leg.set_zorder(110)

    return fig, ax
def relative_vorticity_from_uv(u, v, lon, lat):
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dv_dlon = np.gradient(v, lon, axis=1)
    du_dlat = np.gradient(u, lat, axis=0)
    lat_2d = np.repeat(lat[:, None], len(lon), axis=1)
    R = 6371000
    meters_per_deg_lon = np.deg2rad(1)*R*np.cos(np.deg2rad(lat_2d))
    meters_per_deg_lat = np.deg2rad(1)*R
    meters_per_deg_lon = np.where(np.abs(meters_per_deg_lon) < 1e-12, np.nan, meters_per_deg_lon)
    dvdx = dv_dlon / meters_per_deg_lon
    dudy = du_dlat / meters_per_deg_lat
    return dvdx - dudy

def shear_from_uv(u, v, lon, lat):
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dv_dlon = np.gradient(v, lon, axis=1)
    du_dlat = np.gradient(u, lat, axis=0)
    lat_2d = np.repeat(lat[:, None], len(lon), axis=1)
    R = 6371000
    meters_per_deg_lon = np.deg2rad(1)*R*np.cos(np.deg2rad(lat_2d))
    meters_per_deg_lat = np.deg2rad(1)*R
    meters_per_deg_lon = np.where(np.abs(meters_per_deg_lon) < 1e-12, np.nan, meters_per_deg_lon)
    dvdx = dv_dlon / meters_per_deg_lon
    dudy = du_dlat / meters_per_deg_lat
    return dvdx + dudy

def normal_strain_from_uv(u, v, lon, lat):
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    du_dlon = np.gradient(u, lon, axis=1)
    dv_dlat = np.gradient(v, lat, axis=0)
    lat_2d = np.repeat(lat[:, None], len(lon), axis=1)
    R = 6371000
    meters_per_deg_lon = np.deg2rad(1)*R*np.cos(np.deg2rad(lat_2d))
    meters_per_deg_lat = np.deg2rad(1)*R
    meters_per_deg_lon = np.where(np.abs(meters_per_deg_lon) < 1e-12, np.nan, meters_per_deg_lon)

    dudx = du_dlon / meters_per_deg_lon
    dvdy = dv_dlat / meters_per_deg_lat
    return dudx - dvdy
Omega = 7.292e-5  # Earth's rotation rate (s^-1)

def coriolis_parameter(lat_deg):
    phi = np.deg2rad(lat_deg)
    return 2 * Omega * np.sin(phi)

def OW_from_uv(u, v, lon, lat):
    strain = np.sqrt(normal_strain_from_uv(u, v, lon, lat)**2 + shear_from_uv(u, v, lon, lat)**2)
    vorticity = relative_vorticity_from_uv(u, v, lon, lat)
    OW = strain**2 - vorticity**2

def plot_adcp_satellite_map(
    curr_day,
    sat,
    sat_mag,
    u_adcp,
    v_adcp,
    ugos_track,
    vgos_track,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    longs_bato,
    lat_bato,
):
    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    c = ax.contourf(
        sat["lon"],
        sat["lat"],
        sat_mag,
        levels=np.arange(0, 2.5, 0.1),
        cmap="Spectral_r",
        transform=crs,
    )

    ax.scatter(
        longs_bato[::50],
        lat_bato[::50],
        s=1,
        color="black",
        transform=crs,
        zorder=15,
    )

    step = 50

    # =========================
    # ADCP vectors
    # =========================
    q_adcp = ax.quiver(
        curr_day.lon_adcp.values[::step],
        curr_day.lat_adcp.values[::step],
        u_adcp[::step],
        v_adcp[::step],
        transform=crs,
        scale=3,
        width=0.005,
        color="blueviolet",
        zorder=100,
    )

    # =========================
    # Satellite vectors
    # =========================
    q_sat = ax.quiver(
        curr_day.lon_adcp.values[::step],
        curr_day.lat_adcp.values[::step],
        ugos_track.values[::step],
        vgos_track.values[::step],
        transform=crs,
        scale=3,
        width=0.005,
        color="orange",
        zorder=100,
    )

    # =========================
    # Eddies
    # =========================
    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    # =========================
    # Colorbar
    # =========================
    cb = fig.colorbar(c, ax=ax, shrink=0.6)
    cb.set_label("Geostrophic Velocity [m/s]")

    # =========================
    # Reference arrow (0.5 m/s)
    # =========================
    # =========================
    # Reference arrow (0.5 m/s)
    # =========================
    ax.quiverkey(
        q_adcp,   # or q_sat
        X=1.1,
        Y=0.9,
        U=0.2,
        label="0.2 m/s",
        fontproperties={"size": 8},
        labelpos="E",
        color="black")
        
    # =========================
    # LEGEND (manual handles)
    # =========================
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color="blueviolet", lw=2, label="ADCP"),
        Line2D([0], [0], color="orange", lw=2, label="Satellite"),
    ]

    leg = ax.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize="small",
        frameon=True,
    )
    leg.set_zorder(110) 
    return fig, ax

def get_x_dist_mask(new_adcp, x_dist, dist_min, dist_max):
    """
    Select ADCP data using the cumulative distance along the
    original cruise track.

    Parameters
    ----------
    lon_track, lat_track : 1D arrays
        Original ship track.
    x_dist : 1D array
        Cumulative distance along the original track (m).
    """

    mask = ((x_dist / 1000) >= dist_min) & ((x_dist / 1000) <= dist_max)


    return mask

import numpy as np
from shapely.geometry import LineString, Point, MultiPoint

def get_eddy_adcp_intersections(eddy_contours, eddy_index, adcp, x_dist):
    """
    Find intersections between an eddy contour and an ADCP track.

    Parameters
    ----------
    eddy_contours : GeoDataFrame
        Eddy contours with columns including 'eddy_index' and 'geometry'.
    eddy_index : int
        Eddy index to select.
    lon_adcp, lat_adcp : array-like
        ADCP track longitude and latitude.
    x_dist : array-like
        Cumulative along-track distance (meters) corresponding to ADCP points.

    Returns
    -------
    intersections : list of dict
        Each dictionary contains:
        - lon: longitude of intersection
        - lat: latitude of intersection
        - x_dist: along-track distance (km) at intersection
    """
    lon_adcp = adcp.lon_adcp.values
    lat_adcp = adcp.lat_adcp.values
    # Select eddy contour
    eddy = eddy_contours.loc[
        eddy_contours.eddy_index == eddy_index,
        "geometry"
    ].iloc[0]
    # ADCP track as line
    track = LineString(zip(lon_adcp, lat_adcp))

    # Find intersection
    inter = track.intersection(eddy)

        # Handle different intersection geometries
    if inter.is_empty:
        return []

    if isinstance(inter, Point):
        points = [inter]

    elif isinstance(inter, MultiPoint):
        points = list(inter.geoms)

    elif isinstance(inter, LineString):
        # Track overlaps the contour: use the endpoints of the overlap
        points = [
            Point(inter.coords[0]),
            Point(inter.coords[-1])
        ]

    elif hasattr(inter, "geoms"):
        # GeometryCollection, MultiLineString, etc.
        points = [
            geom for geom in inter.geoms
            if isinstance(geom, Point)
        ]

    else:
        points = []

    results = []

    for p in points:
        lon_int = p.x
        lat_int = p.y

        # Find closest ADCP point
        dist = np.sqrt(
            (np.asarray(lon_adcp) - lon_int)**2 +
            (np.asarray(lat_adcp) - lat_int)**2
        )

        idx = np.argmin(dist)

        results.append({
            "lon": lon_int,
            "lat": lat_int,
            "x_dist": x_dist[idx] / 1000  # km
        })

    return results

import numpy as np
from shapely.geometry import Point, LineString

def determined_x_dist_of_center(eddy_centers, eddy_index, adcp, x_dist):
    """
    Determine the ADCP track location closest to an eddy center.

    Parameters
    ----------
    eddy_centers : GeoDataFrame
        Eddy centers with columns 'eddy_index' and 'geometry' (Point geometries).
    eddy_index : int
        Eddy index to select.
    lon_adcp, lat_adcp : array-like
        ADCP track longitude and latitude.
    x_dist : array-like
        Cumulative along-track distance (meters).

    Returns
    -------
    dict
        Dictionary containing:
        - lon: longitude of closest point on ADCP track
        - lat: latitude of closest point on ADCP track
        - x_dist: along-track distance (km)
    """
    lon_adcp = adcp.lon_adcp.values
    lat_adcp = adcp.lat_adcp.values
    # Select eddy center
    center = eddy_centers.loc[
        eddy_centers.eddy_index == eddy_index,
        "geometry"
    ].iloc[0]

    center_point = Point(center.x, center.y)

    # Create ADCP track
    track = LineString(zip(lon_adcp, lat_adcp))

    # Closest point on track to eddy center
    closest_point = track.interpolate(track.project(center_point))

    lon_closest = closest_point.x
    lat_closest = closest_point.y

    # Compute cumulative distance along track in coordinate space
    segment_lengths = np.sqrt(
        np.diff(lon_adcp)**2 +
        np.diff(lat_adcp)**2
    )

    cum_length = np.concatenate([[0], np.cumsum(segment_lengths)])

    # Distance along track of closest point
    projected_distance = track.project(closest_point)

    # Interpolate real x_dist (meters)
    x_dist_closest = np.interp(
        projected_distance,
        cum_length,
        x_dist
    )

    return {
        "lon": lon_closest,
        "lat": lat_closest,
        "x_dist": x_dist_closest / 1000  # km
    }
import numpy as np

def orthogonal_velocity(adcp, x_dist):
    """
    Calculate velocity component perpendicular to an ADCP track,
    using along-track distance for derivatives.

    Parameters
    ----------
    adcp : xarray.Dataset
        ADCP dataset containing:
        lon_adcp(time), lat_adcp(time),
        uvel_adcp(time, depth), vvel_adcp(time, depth)

    x_dist : array-like
        Cumulative along-track distance in meters.

    Returns
    -------
    v_orth : np.ndarray
        Velocity perpendicular to the track (time, depth).
    """
    
    u = adcp.uvel_adcp.values
    v = adcp.vvel_adcp.values

    lon = adcp.lon_adcp.values
    lat = adcp.lat_adcp.values

    # Convert distance to km (or keep meters)
    x = x_dist / 1000

    # Track derivatives with respect to distance
    dlon_dx = np.gradient(lon, x)
    dlat_dx = np.gradient(lat, x)

    # Normalize tangent vector
    norm = np.sqrt(dlon_dx**2 + dlat_dx**2)

    tx = dlon_dx / norm
    ty = dlat_dx / norm

    # Normal vector (left of track)
    nx = -ty
    ny = tx

    # Expand over depth dimension
    nx = nx[:, None]
    ny = ny[:, None]

    # Project velocity onto normal direction
    v_orth = u * nx + v * ny

    return v_orth

def plot_adt_map_for_adcp(
    sat,
    lon,
    lat,
    anti,
    cyc,
    extent,
    title,
    proj,
    crs,
    ft,
    vmin,
    vmax,
    lon_subset,
    lat_subset,
    depth_subset,
    x_dist,
    dx=50
):

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection=proj)

    setup_map(ax, extent, title, crs, ft)

    bathy = ax.contour(
        lon_subset,
        lat_subset,
        depth_subset,
        colors="grey",
        linewidths=0.8,
        transform=crs,
    )

    ax.clabel(bathy, fontsize=4)

    c = ax.contourf(
        sat["lon"],
        sat["lat"],
        sat["adt"][0],
        levels=np.arange(vmin, vmax, 0.01),
        cmap="RdBu_r",
        transform=crs,
    )
    target_dists = np.arange(np.nanmin(x_dist), np.nanmax(x_dist), dx)  # km

    x_km = x_dist

    for d in target_dists:
        idx = np.argmin(np.abs(x_km - d))

        ax.text(
            lon[idx],
            lat[idx],
            f"{int(d)} km",
            transform=crs,
            fontsize=3,
            color="black",
            ha="left",
            va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="black",
                linewidth=0.2,
                alpha=0.8,
            ),
            zorder=120,
        )

    plot_eddies(ax, anti, "red", crs)
    plot_eddies(ax, cyc, "blue", crs)

    ax.scatter(
        lon,
        lat,
        s=0.5,
        color="green",
        transform=crs,
        zorder=90,
    )


    cb = fig.colorbar(c, ax=ax, shrink=0.6)
    cb.set_label("ADT [m]")

    return fig, ax
