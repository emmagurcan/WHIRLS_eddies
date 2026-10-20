"""
Track one TOEddies eddy day by day through the daily GeoJSON files (WHIRLS, Agulhas / Cape Basin)
and plot its trajectory on the bathymetry, as in the old Labrador Sea (MSM40) script.

How it works
------------
1. You give a starting point (lon, lat), a date range and the eddy type.
2. Each day, among the eddy centers of that type in the TOEddies files, it takes the one
   closest to the previous position (first day: closest to your starting point).
3. If the closest center is further than MAX_JUMP_KM, that day is skipped (the eddy was
   not detected, or the nearest one is a different eddy). After MAX_MISSED_DAYS skipped
   days in a row, tracking stops.
4. The track goes forward in time if START_DATE < END_DATE, backward if START_DATE > END_DATE.

To choose an eddy: set PICK_MODE = True and run. It saves a map of all eddies on START_DATE
with their numbers and prints their centers. Copy the center you want into START_LON/START_LAT,
set PICK_MODE = False and run again.

Outputs (in OUT_DIR)
--------------------
- tracked_eddy_<name>.png : track on bathymetry, centers coloured by date,
                            max-speed contour every CONTOUR_EVERY days
- tracked_eddy_<name>.csv : date, center, distance moved, radius of the max-speed contour
- pick_eddies_<date>.png  : (PICK_MODE) numbered eddies on the start date
"""

import os
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib
if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# =============================================================================
# Settings — choose the eddy here
# =============================================================================
NAME = "eddy_A"              # used in titles and file names
EDDY_TYPE = 1                # 1 = anticyclone, -1 = cyclone
START_LON, START_LAT = 11.6, -37.4
START_DATE = "2026-05-14"    # the day where you know which eddy it is
END_DATE = "2026-07-27"      # earlier than START_DATE -> tracks backward in time

PICK_MODE = True            # True: only draw the numbered eddies of START_DATE to choose one

MAX_JUMP_KM = 60.0           # max distance between consecutive daily centers
MAX_MISSED_DAYS = 3          # stop after this many consecutive days without a match
CONTOUR_EVERY = 3            # draw the max-speed contour every N tracked days (0 = never)

# Map extent: None -> automatic around the track
EXTENT = None                # e.g. [8, 25, -41, -34]
PAD_DEG = 1.5

# =============================================================================
# Paths
# =============================================================================
BASE = Path("/home/geosciences/Documents/Gurcan")
DIR_EDDIES_BETTER = BASE / "data/TOEDDIES/TOEDDIES_better_ADT"   # until 2026-07-22
DIR_EDDIES_STD = BASE / "data/TOEDDIES/TOEDDIES"                  # after that
FNAME_BATHY = BASE / "WHIRLS_eddies/REP/data/bathymetry/etopo2.nc"
OUT_DIR = BASE / "WHIRLS_eddies/figures/eddy_tracks"
OUT_DIR.mkdir(parents=True, exist_ok=True)




# =============================================================================
# Reading the daily TOEddies files
# =============================================================================
def eddy_file(date, kind, part):
    """Path of the <kind>_<part> GeoJSON for this date: better-ADT folder first, then the
    standard one. kind = 'anticyclones'|'cyclones', part = 'centers'|'maxspeed'|'outcontour'."""
    d = pd.Timestamp(date)
    for folder, tag in ((DIR_EDDIES_BETTER, d.strftime("%Y%m%d")),
                        (DIR_EDDIES_STD, d.strftime("%Y-%m-%d")),
                        (DIR_EDDIES_STD, d.strftime("%Y%m%d"))):
        for f in sorted(glob.glob(str(folder / f"*{tag}*"))):
            name = Path(f).name.lower()
            is_anti = "anticyclon" in name
            # careful: "cyclones" is also inside "anticyclones"
            if (kind == "anticyclones") != is_anti:
                continue
            if part in name and name.endswith((".geojson", ".json", ".shp")):
                return f, folder.name
    return None, None


def load_day(date, kind):
    """Centers (lon, lat) and max-speed contours for one day and one eddy type."""
    fc, src = eddy_file(date, kind, "centers")
    if fc is None:
        return None
    cen = gpd.read_file(fc)
    fm, _ = eddy_file(date, kind, "maxspeed")
    ms = gpd.read_file(fm) if fm else None
    ok = cen.geometry.notna() & ~cen.geometry.is_empty
    lon = np.full(len(cen), np.nan)
    lat = np.full(len(cen), np.nan)
    lon[ok.values] = cen.geometry[ok].x.values
    lat[ok.values] = cen.geometry[ok].y.values
    lon = np.where(lon > 180, lon - 360, lon)
    return {"lon": lon, "lat": lat, "centers": cen, "maxspeed": ms, "source": src}


def contour_xy(geom):
    """Outline vertices of a (Multi)Polygon / LineString, NaN vertices removed."""
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type == "Polygon":
        xy = np.asarray(geom.exterior.coords)
    elif geom.geom_type == "MultiPolygon":
        xy = np.asarray(max(geom.geoms, key=lambda g: g.area).exterior.coords)
    elif geom.geom_type in ("LineString", "LinearRing"):
        xy = np.asarray(geom.coords)
    elif geom.geom_type == "MultiLineString":
        xy = np.asarray(max(geom.geoms, key=lambda g: g.length).coords)
    else:
        return None
    xy = xy[:, :2][np.isfinite(xy[:, :2]).all(axis=1)]
    xy[:, 0] = np.where(xy[:, 0] > 180, xy[:, 0] - 360, xy[:, 0])
    return xy if len(xy) > 2 else None


def maxspeed_contour(day, i):
    """Max-speed contour of eddy i. The TOEddies files list eddies in the same order in the
    centers and maxspeed files; if not, use the contour that contains the center."""
    ms = day["maxspeed"]
    if ms is None or len(ms) == 0:
        return None
    if len(ms) == len(day["centers"]):
        xy = contour_xy(ms.geometry.iloc[i])
        if xy is not None:
            return xy
    p = day["centers"].geometry.iloc[i]
    for g in ms.geometry:
        if g is not None and not g.is_empty and g.geom_type.endswith("Polygon") and g.contains(p):
            return contour_xy(g)
    return None


# =============================================================================
# Geometry helpers
# =============================================================================
def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, (lon1, lat1, lon2, lat2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(a))


def equivalent_radius_km(xy):
    """Radius of the circle with the same area as the contour (local flat projection)."""
    if xy is None:
        return np.nan
    lat0 = np.deg2rad(np.mean(xy[:, 1]))
    x = np.deg2rad(xy[:, 0]) * R_EARTH_KM * np.cos(lat0)
    y = np.deg2rad(xy[:, 1]) * R_EARTH_KM
    area = 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    return np.sqrt(area / np.pi)


# =============================================================================
# Map background (same style as the MSM40 figure)
# =============================================================================
def setup_map(fig, extent, sub=111):
    proj = ccrs.Mercator(central_longitude=(extent[0] + extent[1]) / 2)
    ax = fig.add_subplot(sub, projection=proj)
    ax.set_extent(extent, crs=crs)
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="gray", zorder=3)
    ax.coastlines(resolution="10m", color="k", linewidth=0.6, zorder=4)
    gl = ax.gridlines(linewidth=1, color="black", alpha=0.2, linestyle="--",
                      draw_labels=True, zorder=6)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": ft}
    b = xr.open_dataset(FNAME_BATHY)
    x0, x1, y0, y1 = extent
    topo = b["topo"].sel(lon=slice(x0 - 1, x1 + 1), lat=slice(y0 - 1, y1 + 1))
    c = ax.contourf(topo.lon, topo.lat, topo, np.arange(-5500, 1, 500),
                    transform=crs, cmap="Blues_r", zorder=0, alpha=0.85, extend="min")
    return ax, c


# =============================================================================
# Pick mode: show the numbered eddies of the start date
# =============================================================================
def pick_eddies():
    day = load_day(START_DATE, KIND)
    if day is None:
        raise FileNotFoundError(f"No {KIND} centers file for {START_DATE}")
    extent = EXTENT or [START_LON - 6, START_LON + 6, START_LAT - 4, START_LAT + 4]
    x0, x1, y0, y1 = extent
    inside = np.where((day["lon"] >= x0) & (day["lon"] <= x1)
                      & (day["lat"] >= y0) & (day["lat"] <= y1))[0]
    print(f"{KIND} on {START_DATE} ({day['source']}), inside {extent}:")
    print(f"{'#':>4} {'lon':>9} {'lat':>9} {'dist to start (km)':>20}")
    for i in inside:
        dkm = haversine_km(START_LON, START_LAT, day["lon"][i], day["lat"][i])
        print(f"{i:>4} {day['lon'][i]:9.3f} {day['lat'][i]:9.3f} {dkm:20.0f}")

    fig = plt.figure(figsize=(9, 7), dpi=150)
    ax, _ = setup_map(fig, extent)
    col = "r" if EDDY_TYPE == 1 else "b"
    for i in inside:
        xy = maxspeed_contour(day, i)
        if xy is not None:
            ax.plot(xy[:, 0], xy[:, 1], color=col, lw=1.5, transform=crs, zorder=10)
        ax.scatter(day["lon"][i], day["lat"][i], s=12, c="k", transform=crs, zorder=11)
        ax.text(day["lon"][i] + 0.08, day["lat"][i] + 0.08, str(i), fontsize=8,
                transform=crs, zorder=12, weight="bold")
    ax.scatter(START_LON, START_LAT, marker="*", s=150, c="yellow", edgecolors="k",
               transform=crs, zorder=13, label="START_LON/LAT")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"{KIND.capitalize()} on {START_DATE}: pick one, copy its center", fontsize=11)
    out = OUT_DIR / f"pick_eddies_{pd.Timestamp(START_DATE):%Y%m%d}_{KIND}.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# =============================================================================
# Tracking
# =============================================================================
def track():
    t0, t1 = pd.Timestamp(START_DATE), pd.Timestamp(END_DATE)
    step = pd.Timedelta(days=1) if t1 >= t0 else pd.Timedelta(days=-1)
    dates = pd.date_range(t0, t1, freq=step)

    cur_x, cur_y = START_LON, START_LAT
    rows, missed = [], 0
    for date in dates:
        day = load_day(date, KIND)
        if day is None:
            print(f"{date:%Y-%m-%d}: no {KIND} file, skipped")
            missed += 1
        else:
            d = haversine_km(cur_x, cur_y, day["lon"], day["lat"])
            if np.all(np.isnan(d)):
                missed += 1
                continue
            i = int(np.nanargmin(d))
            limit = MAX_JUMP_KM * (missed + 1)          # allow a longer jump after gaps
            if d[i] > limit:
                print(f"{date:%Y-%m-%d}: nearest {KIND[:-1]} is {d[i]:.0f} km away "
                      f"(> {limit:.0f} km), skipped")
                missed += 1
            else:
                xy = maxspeed_contour(day, i)
                rows.append(dict(date=date, lon=day["lon"][i], lat=day["lat"][i],
                                 step_km=d[i] if rows else np.nan,
                                 r_maxspeed_km=equivalent_radius_km(xy),
                                 index_in_file=i, source=day["source"], contour=xy))
                print(f"{date:%Y-%m-%d}: #{i:<4} {day['lon'][i]:8.3f} {day['lat'][i]:8.3f}"
                      f"  moved {d[i]:5.1f} km  Rmax {rows[-1]['r_maxspeed_km']:5.1f} km"
                      f"  ({day['source']})")
                cur_x, cur_y = day["lon"][i], day["lat"][i]
                missed = 0
        if missed > MAX_MISSED_DAYS:
            print(f"Lost the eddy after {MAX_MISSED_DAYS} days without a match: stop.")
            break

    tr = pd.DataFrame(rows)
    if tr.empty:
        raise RuntimeError("No eddy found: check START_LON/LAT, START_DATE and EDDY_TYPE "
                           "(run with PICK_MODE = True).")
    tr = tr.sort_values("date").reset_index(drop=True)
    tr.drop(columns="contour").to_csv(OUT_DIR / f"tracked_eddy_{NAME}.csv", index=False)
    return tr


def plot_track(tr):
    if EXTENT is not None:
        extent = EXTENT
    else:
        xs = np.concatenate([tr.lon.values] + [c[:, 0] for c in tr.contour if c is not None])
        ys = np.concatenate([tr.lat.values] + [c[:, 1] for c in tr.contour if c is not None])
        extent = [xs.min() - PAD_DEG, xs.max() + PAD_DEG, ys.min() - PAD_DEG, ys.max() + PAD_DEG]

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax, c = setup_map(fig, extent)

    tnum = mdates.date2num(tr.date)
    norm = Normalize(vmin=tnum.min(), vmax=tnum.max() if tnum.max() > tnum.min() else tnum.min() + 1)
    cmap = plt.get_cmap("autumn")

    ax.plot(tr.lon, tr.lat, "-", color="k", lw=0.8, transform=crs, zorder=9)
    if CONTOUR_EVERY:
        idx = list(range(0, len(tr), CONTOUR_EVERY))
        if idx[-1] != len(tr) - 1:
            idx.append(len(tr) - 1)
        for k in idx:
            xy = tr.contour[k]
            if xy is not None:
                ax.plot(xy[:, 0], xy[:, 1], color=cmap(norm(tnum[k])), lw=1.2,
                        transform=crs, zorder=8)
    ax.scatter(tr.lon, tr.lat, c=tnum, cmap=cmap, norm=norm, s=14, edgecolors="k",
               linewidths=0.3, transform=crs, zorder=10)
    ax.scatter(START_LON, START_LAT, marker="*", s=120, c="yellow", edgecolors="k",
               transform=crs, zorder=11)

    typ = "Anticyclone" if EDDY_TYPE == 1 else "Cyclone"
    ax.set_title(f"{typ} {NAME}: tracked centers "
                 f"({tr.date.min():%b %d} – {tr.date.max():%b %d, %Y})", fontsize=12)
    ax.legend(handles=[
        mlines.Line2D([], [], marker="o", ls="-", color="k", mfc="orange", ms=5,
                      label="Daily center"),
        mlines.Line2D([], [], color="orange", lw=1.2,
                      label=f"Max-speed contour (every {CONTOUR_EVERY} d)"),
        mlines.Line2D([], [], marker="*", ls="", mfc="yellow", mec="k", ms=10,
                      label="Starting point"),
    ], loc="upper right", fontsize=7, framealpha=0.9).set_zorder(20)

    # depth colorbar under the map, date colorbar on the side (as in the MSM40 figure)
    cax = inset_axes(ax, width="80%", height="3%", loc="lower center",
                     bbox_to_anchor=(0, -0.12, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
    cb = fig.colorbar(c, cax=cax, orientation="horizontal")
    cb.set_label("Depth [m]", fontsize=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb2 = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.6, pad=0.03)
    cb2.ax.yaxis.set_major_locator(mdates.AutoDateLocator(maxticks=6))
    cb2.ax.yaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    cb2.set_label("Date", fontsize=9)

    out = OUT_DIR / f"tracked_eddy_{NAME}.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# =============================================================================
if __name__ == "__main__":
    if PICK_MODE:
        pick_eddies()
    else:
        tr = track()
        print(f"\n{len(tr)} days tracked, {tr.step_km.sum():.0f} km in total, "
              f"mean drift {tr.step_km.mean():.1f} km/day, "
              f"mean Rmax {tr.r_maxspeed_km.mean():.0f} km")
        plot_track(tr)