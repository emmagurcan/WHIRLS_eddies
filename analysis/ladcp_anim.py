"""
Daily animation of ADT + TOEddies contours + LADCP stations over the cruise.

For each day between CRUISE_START and CRUISE_END:
  - ADT of that day (filled contours)
  - anticyclones (red) / cyclones (blue) detected that day
  - LADCP stations done THAT day: large yellow markers with station numbers
  - LADCP stations done on EARLIER days: small faded grey markers (cruise history)
  - stations not yet done are not shown

Frames are saved as PNG in FRAME_DIR, then combined into a GIF and an MP4.
"""

import os
import warnings
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")                     # no window: we only save files
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from ladcp_tools import *                  # load_all_stations, ...
from tools import *                        # load_satellite, load_eddies, plot_eddies, setup_map

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
lon_min, lon_max, lat_min, lat_max = 8, 20, -41, -34
extent = [lon_min, lon_max, lat_min, lat_max]

CRUISE_START = "2026-06-28"
CRUISE_END   = "2026-07-27"

ft = 10
vmin, vmax = -1.5, 1.5
levels = np.arange(vmin, vmax + 0.01, 0.01)

dir_ladcp  = "/home/geosciences/Documents/Gurcan/data/LADCP/traitement_avec_SADCP/"
dir_eddies_better = "/home/geosciences/Documents/Gurcan/data/TOEDDIES/TOEDDIES_better_ADT/"  # until 2026-07-22
dir_eddies_std    = "/home/geosciences/Documents/Gurcan/data/TOEDDIES/TOEDDIES/"             # after that
dir_ssh    = "/home/geosciences/Documents/Gurcan/data/SSH/"
fname_bathy = "/home/geosciences/Documents/Gurcan/WHIRLS_eddies/REP/data/bathymetry/etopo2.nc"

FRAME_DIR = "frames_ladcp_animation"
OUT_GIF   = "ladcp_eddies_animation.gif"
OUT_MP4   = "ladcp_eddies_animation.mp4"
FPS = 2                                   # frames (days) per second

os.makedirs(FRAME_DIR, exist_ok=True)

crs = ccrs.PlateCarree()
proj = ccrs.Mercator(central_longitude=(lon_min + lon_max) / 2)

# ----------------------------------------------------------------------------
# Things loaded ONCE (they don't change from day to day)
# ----------------------------------------------------------------------------
# Bathymetry
bathymetry_data = xr.open_dataset(fname_bathy)
lon_subset = bathymetry_data["lon"].sel(lon=slice(lon_min, lon_max))
lat_subset = bathymetry_data["lat"].sel(lat=slice(lat_min, lat_max))
depth_subset = bathymetry_data["topo"].sel(lon=slice(lon_min, lon_max),
                                           lat=slice(lat_min, lat_max))

# File lists
daily_ssh_files = sorted(os.listdir(dir_ssh))
eddy_files_better = sorted(os.listdir(dir_eddies_better))
eddy_files_std = sorted(os.listdir(dir_eddies_std))


def get_eddy_files(date_compact, date):
    """Use the better-ADT TOEddies when available for that day, else the standard ones."""
    files = [dir_eddies_better + f for f in eddy_files_better if date_compact in f]
    if files:
        return files, "TOEddies (better ADT)"
    files = [dir_eddies_std + f for f in eddy_files_std if date in f]
    if files:
        return files, "TOEddies (standard)"
    return [], None

# LADCP stations
st = load_all_stations(dir_ladcp)
st = st[(st.lon >= lon_min) & (st.lon <= lon_max) &
        (st.lat >= lat_min) & (st.lat <= lat_max)].copy()


def find_time_column(df):
    """Return the name of the column holding the cast time."""
    for cand in ["time", "date", "datetime", "start_time", "time_start", "t"]:
        if cand in df.columns:
            return cand
    for c in df.columns:                       # any datetime-typed column
        if np.issubdtype(df[c].dtype, np.datetime64):
            return c
    raise KeyError(
        "No time column found in the LADCP station table. "
        f"Columns are: {list(df.columns)}. "
        "Add the cast time in load_all_stations() or set TIME_COL by hand."
    )


TIME_COL = find_time_column(st)
st["day"] = pd.to_datetime(st[TIME_COL]).dt.normalize()


# ----------------------------------------------------------------------------
# One frame
# ----------------------------------------------------------------------------
def make_frame(date, frame_idx):
    date_str = date.strftime("%Y-%m-%d")
    date_compact = date.strftime("%Y%m%d")

    # --- ADT of the day ---
    ssh_match = [f for f in daily_ssh_files if date_str in f or date_compact in f]
    if not ssh_match:
        warnings.warn(f"No ADT file for {date_str}, skipping this day.")
        return None
    ssh = load_satellite(dir_ssh + ssh_match[0])

    # --- Eddies of the day ---
    eddy_files, eddy_source = get_eddy_files(date_compact, date_str)

    # Fixed figure size and no bbox_inches="tight" so all frames are identical in size
    fig = plt.figure(figsize=(7, 4.8), dpi=150)
    ax = fig.add_subplot(111, projection=proj)
    setup_map(ax, extent, f"ADT, eddies and LADCP stations: {date_str}", crs, ft)

    bathy = ax.contour(lon_subset, lat_subset, depth_subset,
                       colors="grey", linewidths=0.6, transform=crs)
    ax.clabel(bathy, fontsize=4)

    c = ax.contourf(ssh["lon"], ssh["lat"], ssh["adt"][0],
                    levels=levels, cmap="RdBu_r", extend="both",
                    transform=crs)

    if eddy_files:
        anti, cyc = load_eddies(eddy_files)
        plot_eddies(ax, anti, "red", crs)
        plot_eddies(ax, cyc, "blue", crs)
    else:
        warnings.warn(f"No eddy file for {date_str} in either TOEddies folder.")

    # Small label so you can see which eddy product each frame uses
    if eddy_source:
        ax.text(0.01, 0.01, eddy_source, transform=ax.transAxes, fontsize=5,
                ha="left", va="bottom", zorder=102,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1))

    # --- LADCP stations ---
    past = st[st.day < date]
    today = st[st.day == date]

    ax.scatter(past.lon, past.lat, s=6, color="green", edgecolors="green",
               linewidths=0.3, alpha=0.8, zorder=99, transform=crs)

    ax.scatter(today.lon, today.lat, s=22, color="yellow", edgecolors="black",
               linewidths=0.6, zorder=100, transform=crs)

    for _, s in today.iterrows():
        ax.annotate(str(int(s.station)), xy=(s.lon, s.lat),
                    xycoords=crs._as_mpl_transform(ax),
                    xytext=(3, 3), textcoords="offset points",
                    fontsize=5, zorder=101,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="black", linewidth=0.3, alpha=0.85))

    # --- Legend (proxy handles, so it is the same on every frame) ---
    handles = [
        mlines.Line2D([], [], marker="o", ls="", color="yellow", mec="black",
                      ms=5, label=f"LADCP today ({len(today)})"),
        mlines.Line2D([], [], marker="o", ls="", color="green", mec="green",
                      ms=3, label="LADCP earlier"),
        mlines.Line2D([], [], color="red", label="Anticyclones"),
        mlines.Line2D([], [], color="blue", label="Cyclones"),
    ]
    ax.legend(handles=handles, fontsize="x-small", loc="upper right",
              framealpha=0.85)

    cb = fig.colorbar(c, ax=ax, shrink=0.6, ticks=np.arange(vmin, vmax + 0.01, 0.3))
    cb.set_label("ADT [m]")

    fname = os.path.join(FRAME_DIR, f"frame_{frame_idx:03d}_{date_compact}.png")
    fig.savefig(fname)
    plt.close(fig)
    return fname


# ----------------------------------------------------------------------------
# Loop over the cruise
# ----------------------------------------------------------------------------
dates = pd.date_range(CRUISE_START, CRUISE_END, freq="D")
frames = []
for i, d in enumerate(dates):
    print(f"Frame {i + 1}/{len(dates)}: {d.date()}")
    f = make_frame(d, i)
    if f is not None:
        frames.append(f)

# ----------------------------------------------------------------------------
# Assemble the animation
# ----------------------------------------------------------------------------
import imageio.v2 as imageio

images = [imageio.imread(f) for f in frames]
imageio.mimsave(OUT_GIF, images, duration=1000 / FPS, loop=0)   # duration in ms
print(f"Saved {OUT_GIF}")

try:
    imageio.mimsave(OUT_MP4, images, fps=FPS, macro_block_size=1)
    print(f"Saved {OUT_MP4}")
except Exception as e:
    print(f"MP4 not written ({e}). Install imageio-ffmpeg: pip install imageio-ffmpeg")