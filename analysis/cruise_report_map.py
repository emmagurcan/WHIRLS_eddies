"""
WHIRLS cruise report: overview map (two panels)
==============================================
Both panels: SST (CLS 7-day composite) + TOEddies contours (better ADT) of DATE.

    (a) left : Argo/profiling floats, Seagliders + SeaExplorer, Wave Gliders, XSPAR,
               CTD stations, MVP transects
    (b) right: surface drifter trajectories coloured by speed, deployment positions

Platform positions cover the whole cruise; SST and eddies are a single-day snapshot.
Both maps have the same extent and share one SST colour bar.

Every layer can be switched off in SHOW. Missing files are skipped with a message.
Column names in the CSVs are found automatically (lon/longitude/LON (°E), time/date, ...).
"""
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd

try:
    import cmocean.cm as cmo
    SST_CMAP, SPEED_CMAP = cmo.thermal, cmo.speed
except ImportError:
    SST_CMAP, SPEED_CMAP = "magma", "YlGn"

# =============================================================================
# Paths
# =============================================================================
BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
DIR_SADCP = BASE / "data/SADCP"
SADCP_FILES = {"OS150": DIR_SADCP / "WHIRLS_NB150_TOT_osite.nc",
               "OS38": DIR_SADCP / "WHIRLS_BB38_TOT_osite.nc"}
DIR_SSH = BASE / "data/SSH"
DIR_EDDIES = {"TOEDDIES_better_ADT": BASE / "data/TOEDDIES/TOEDDIES_better_ADT",
              "TOEDDIES": BASE / "data/TOEDDIES/TOEDDIES"}
DIR_TRACKS = BASE / "eddies"
OUT_DIR = BASE / "figures"

DIR_FLOATS = BASE / "data/FLOATS"
DIR_DRIFTERS = BASE / "data/DRIFTERS/drifters"                    # folder of CSV snapshots

DIR_MVP = BASE / "data/MVP_Data_processed"
DIR_SG = BASE / "data/SEAGLIDERS"
DIR_WG = BASE / "data/WAVEGLIDERS"
DIR_XSP = BASE / "data/XSPAR"
DIR_SST = BASE / "data/SST"
DIR_CTD_META = BASE / "data/CTD/ctd_metadata.nc"
FNAME_BATHY = BASE / "data/bathymetry/etopo2.nc"
crs = ccrs.PlateCarree()
DIR_CODE = BASE / "analysis"            # where tools_eg.py lives
sys.path.insert(0, str(DIR_CODE))

# =============================================================================
# CONFIG
# =============================================================================
DATE = "2026-07-13"                                   # SST + eddy snapshot (both panels)
SST_FILE = DIR_SST / f"{DATE.replace('-', '')}_CLS_SST_7d.nc"
EDDY_PRODUCT = "TOEDDIES_better_ADT"                  # falls back to "TOEDDIES" if missing
CRUISE_END = "2026-07-27"

EXTENT = None            # [lon_min, lon_max, lat_min, lat_max]; None = fit all platforms
EXTENT_PAD = 0.6         # deg added around the platforms when EXTENT is None
SST_RANGE = None         # (vmin, vmax) in °C; None = 2nd-98th percentile in the map
DRIFTER_PANEL_SST = "colour"   # "colour" or "grey": grey SST makes the speed colours pop

SHOW = dict(sst=True, eddies=True, eddy_outer=True, bathy=True,
            floats=True, gliders=True, wavegliders=True, xspar=True, ctd=True, mvp=True,
            drifters=True)
LABEL_MVP = True         # write the transect number at the start of each MVP line
BATHY_LEVELS = [-3000, -1000]
MAX_GAP_KM = 30.0        # break a track where two fixes are further apart (recovery/redeploy)

# Drifters
FIRST_LAUNCH = "2026-07-01 05:00"  # fixes before this are bench tests
DRIFTER_END = CRUISE_END           # last fix drawn (None = whole record)
MAX_DT_H = 6.0                     # no segment across a gap longer than this
MAX_SPEED = 3.0                    # m/s: faster segments are GPS spikes -> dropped
SPEED_SMOOTH_N = 5                 # running mean over this many segments (1 = raw)
SPEED_LIM = None                   # (vmin, vmax) m/s; None = 0 .. 98th percentile

# Autonomous platforms
GLIDER_GLOB = "sg*_track.csv"
SEAEXPLORER_FILE = DIR_SG / "seaexplorer.csv"
WAVEGLIDERS = [(DIR_WG / "melktert_track.csv", "Wave Glider Melktert"),
               (DIR_WG / "vetkoek_track.csv", "Wave Glider Vetkoek")]
XSPAR_FILE = DIR_XSP / "xspar_xeos_track.csv"
FLOAT_PATTERNS = ("*.csv", "*.nc", "*.txt")    # files in DIR_FLOATS

# Colours (colour-blind safe). Gliders = blues/greens, wave gliders = oranges.
GLIDER_COLOURS = ["#56B4E9", "#0072B2", "#009E73", "#00A0A0", "#7FD1C2"]
WAVEGLIDER_COLOURS = ["#E69F00", "#D55E00"]
XSPAR_COLOUR = "#CC79A7"
FLOAT_COLOUR = "#882255"
MVP_COLOUR = "black"
CTD_FACE = "white"
DRIFTER_FACE = "#F0E442"
AC_COLOUR, C_COLOUR = "black", "white"      # anticyclones in warm water, cyclones in cold

LAND_RES = "50m"         # "10m" is sharper but much slower (and downloads on first use)
SST_COARSEN = None       # e.g. 2 or 3 to average SST cells if the file is very high-res
MVP_CACHE = OUT_DIR / "mvp_cast_positions.csv"   # delete it to re-read the MVP CSVs
MVP_SKIP = ("essais",)   # MVP folders left out of the map (test casts)
TYPE_SADCP = "OS38"      # SADCP used to fill missing MVP positions (as in mvp.ipynb)

FIGSIZE = (17, 9)
FONTSIZE = 10
OUT_NAME = f"WHIRLS_overview_map_2panel_{DATE.replace('-', '')}"

plt.rcParams.update({"font.size": FONTSIZE, "axes.titlesize": FONTSIZE + 2})


# =============================================================================
# Helpers
# =============================================================================
def halo(lw, colour="white"):
    return [pe.Stroke(linewidth=lw, foreground=colour), pe.Normal()]


def find_col(df, exact, contains=None, required=True):
    """Column whose name matches one of `exact` (case-insensitive), else contains a key."""
    low = {str(c).lower().strip(): c for c in df.columns}
    for k in exact:
        if k.lower() in low:
            return low[k.lower()]
    for k in contains or []:
        for cl, c in low.items():
            if k in cl:
                return c
    if required:
        raise KeyError(f"No column like {exact} in {list(df.columns)}")
    return None


def lonlat_cols(df):
    lon = find_col(df, ["lon", "longitude", "LON (°E)", "long", "lon_deg"], ["lon"])
    lat = find_col(df, ["lat", "latitude", "LAT (°N)", "lat_deg"], ["lat"])
    return lon, lat


def time_col(df):
    return find_col(df, ["time", "datetime", "date", "timestamp", "TIME", "time_utc", "juld"],
                    ["time", "date"], required=False)


def id_col(df):
    return find_col(df, ["id", "drifter_id", "float_id", "platform_number", "wmo", "drifter",
                         "imei", "name", "platform", "buoy", "device", "deviceid", "serial"],
                    ["imei", "wmo"], required=False)


def clean_lonlat(lon, lat):
    lon = pd.to_numeric(pd.Series(np.ravel(lon)), errors="coerce").values.astype(float)
    lat = pd.to_numeric(pd.Series(np.ravel(lat)), errors="coerce").values.astype(float)
    lon = np.where(lon > 180, lon - 360, lon)
    bad = (np.abs(lat) > 90) | (np.abs(lon) > 180) | ((lon == 0) & (lat == 0))
    lon[bad], lat[bad] = np.nan, np.nan
    return lon, lat


def break_gaps(lon, lat, max_km=MAX_GAP_KM):
    """Insert NaN where consecutive fixes are > max_km apart, so no straight jumps are drawn."""
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]
    if lon.size < 2:
        return lon, lat
    dx = np.deg2rad(np.diff(lon)) * np.cos(np.deg2rad(0.5 * (lat[1:] + lat[:-1])))
    d = 6371.0 * np.hypot(dx, np.deg2rad(np.diff(lat)))
    cut = np.flatnonzero(d > max_km) + 1
    return np.insert(lon, cut, np.nan), np.insert(lat, cut, np.nan)


def read_track(path):
    """Platform track CSV -> (lon, lat), time-sorted, gaps broken. None if missing."""
    path = Path(path)
    if not path.exists():
        print(f"  ! missing {path}")
        return None
    df = pd.read_csv(path)
    clon, clat = lonlat_cols(df)
    tc = time_col(df)
    if tc is not None:
        df = df.assign(_t=pd.to_datetime(df[tc], errors="coerce")).sort_values("_t")
    lon, lat = clean_lonlat(df[clon], df[clat])
    print(f"  {path.name}: {np.isfinite(lon).sum()} fixes")
    return break_gaps(lon, lat)


_T0 = [time.time()]


def tick(msg):
    """Print msg with the time since the previous tick (to see which step is slow)."""
    now = time.time()
    print(f"  [{now - _T0[0]:6.1f} s] {msg}", flush=True)
    _T0[0] = now


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]


# =============================================================================
# Loaders
# =============================================================================
def load_sst(path, date):
    ds = xr.open_dataset(path)
    names = ["analysed_sst", "sst", "SST", "sea_surface_temperature",
             "adjusted_sea_surface_temperature", "thetao"]
    var = next((n for n in names if n in ds.data_vars), None)
    if var is None:
        var = next(v for v in ds.data_vars if ds[v].ndim >= 2)
    da = ds[var]
    lon = next(c for c in ["lon", "longitude", "LONGITUDE", "x"] if c in da.coords or c in ds)
    lat = next(c for c in ["lat", "latitude", "LATITUDE", "y"] if c in da.coords or c in ds)
    if "time" in da.dims:
        da = da.sel(time=np.datetime64(date), method="nearest")
    da = da.squeeze(drop=True)
    if float(da[lon].max()) > 180:
        da = da.assign_coords({lon: ((da[lon] + 180) % 360) - 180})
    da = da.sortby(lon).sortby(lat)
    units = str(da.attrs.get("units", "")).lower()
    if SST_COARSEN:
        da = da.coarsen({lon: SST_COARSEN, lat: SST_COARSEN}, boundary="trim").mean()
    print(f"  SST grid: {da.sizes[lat]} x {da.sizes[lon]}")
    if units in ("k", "kelvin", "degk") or float(da.mean()) > 100:
        da = da - 273.15
    print(f"  SST: '{var}' from {Path(path).name}")
    return da.rename({lon: "lon", lat: "lat"})


def get_eddy_files(date, product=EDDY_PRODUCT):
    """Shapefiles of that day; better-ADT first, standard TOEddies as fallback."""
    tags = (date.replace("-", ""), date)
    for key in [product] + [k for k in DIR_EDDIES if k != product]:
        d = DIR_EDDIES[key]
        if not d.exists():
            continue
        files = [f for f in d.iterdir() if any(t in f.name for t in tags)]
        files = [f for f in files if f.suffix == ".shp"] or files
        if files:
            if key != product:
                print(f"  ! no {product} files for {date}, using {key}")
            return sorted(files), key
    print(f"  ! no TOEddies files for {date}")
    return [], None


def load_eddies(files):
    out = {"anti": {}, "cyc": {}}
    for f in files:
        kind = "anti" if "anti" in f.name.lower() else "cyc"
        for part in ("maxspeed", "outcontour"):
            if part in f.name.lower() and part not in out[kind]:
                out[kind][part] = gpd.read_file(f)
    return out


def load_ctd(path):
    ds = xr.open_dataset(path)
    lon = next(n for n in ["longitude", "lon", "lon_start", "LONGITUDE"] if n in ds)
    lat = next(n for n in ["latitude", "lat", "lat_start", "LATITUDE"] if n in ds)
    lo, la = clean_lonlat(ds[lon].values, ds[lat].values)
    print(f"  CTD: {np.isfinite(lo).sum()} casts")
    return lo, la


# ---- MVP: same steps as mvp.ipynb, using the functions in tools_eg.py ----
_SADCP = {}


def get_sadcp():
    """SADCP loaded once, exactly as in mvp.ipynb (cell 5)."""
    if "ds" not in _SADCP:
        from tools_eg import cascade_add_time
        sadcp = cascade_add_time(xr.open_dataset(SADCP_FILES[TYPE_SADCP]))
        sadcp = sadcp.sortby("time")
        _SADCP["ds"] = sadcp.isel(time=~sadcp.indexes["time"].duplicated())
    return _SADCP["ds"]


def load_mvp(dir_mvp):
    """{transect name: (lon, lat)}, one point per cast.

    Per transect, as in mvp.ipynb:
      1. read every CSV of the folder, CAST = file index      (cell 4)
      2. fill NaN LAT/LON from the SADCP navigation          (cell 5)
      3. split_by_cast_id + fill_cast_positions (tools_eg)    (process_transect)
    The cast position is the mean LAT/LON of the cast (as build_section).
    """
    from tools_eg import (COLS, sadcp_for_transect, sadcp_on_mvp_times,
                          split_by_cast_id, fill_cast_positions)
    if MVP_CACHE.exists():
        c = pd.read_csv(MVP_CACHE)
        print(f"  MVP: positions from cache {MVP_CACHE.name}")
        return {k: (g["lon"].values, g["lat"].values)
                for k, g in sorted(c.groupby("transect", sort=False), key=lambda kg: natural_key(kg[0]))}

    transects = {}
    folders = sorted([d for d in Path(dir_mvp).iterdir() if d.is_dir()
                      and not any(k.lower() in d.name.lower() for k in MVP_SKIP)],
                     key=lambda d: natural_key(d.name))
    for transect_path in folders:
        # ---- cell 4 ----
        files = sorted(Path(transect_path).glob("*.csv"))
        if not files:
            continue
        dfs = []
        for k, f in enumerate(files):
            df = pd.read_csv(f)
            df["CAST"] = k          # one processed CSV file = one MVP cast
            dfs.append(df)
        mvp_trans = pd.concat(dfs, ignore_index=True)

        # ---- cell 5 (only when some positions are missing; lat/lon only) ----
        if mvp_trans[COLS["lat"]].isna().any() or mvp_trans[COLS["lon"]].isna().any():
            sadcp_tr = sadcp_for_transect(mvp_trans, get_sadcp())
            sadcp_mvp = sadcp_on_mvp_times(mvp_trans, sadcp_tr, max_gap="30min",
                                           keep=("LATITUDE", "LONGITUDE"))
            mvp_trans[COLS["lon"]] = mvp_trans[COLS["lon"]].fillna(
                pd.Series(sadcp_mvp["LONGITUDE"].values, index=mvp_trans.index))
            mvp_trans[COLS["lat"]] = mvp_trans[COLS["lat"]].fillna(
                pd.Series(sadcp_mvp["LATITUDE"].values, index=mvp_trans.index))

        # ---- process_transect: casts + their positions ----
        try:
            casts = fill_cast_positions(split_by_cast_id(mvp_trans, cast_col="CAST"))
        except ValueError as e:
            print(f"  ! MVP {transect_path.name} skipped ({e})")
            continue
        lon = np.array([c[COLS["lon"]].mean() for c in casts])
        lat = np.array([c[COLS["lat"]].mean() for c in casts])
        lon = np.where(lon > 180, lon - 360, lon)
        transects[transect_path.name] = (lon, lat)
        tick(f"MVP {transect_path.name}: {len(casts)} casts")

    print(f"  MVP: {len(transects)} transects")
    if transects:
        pd.concat([pd.DataFrame({"transect": k, "lon": v[0], "lat": v[1]})
                   for k, v in transects.items()]).to_csv(MVP_CACHE, index=False)
    return transects


def _read_any(f):
    """DataFrame with lon, lat, time (+ id) from a CSV/TXT or an Argo-style NetCDF."""
    if f.suffix == ".nc":
        ds = xr.open_dataset(f)
        lon = next(n for n in ["LONGITUDE", "longitude", "lon"] if n in ds)
        lat = next(n for n in ["LATITUDE", "latitude", "lat"] if n in ds)
        tn = next((n for n in ["JULD", "TIME", "time", "JULD_LOCATION"] if n in ds), None)
        df = pd.DataFrame({"lon": np.ravel(ds[lon].values), "lat": np.ravel(ds[lat].values)})
        if tn is not None and ds[tn].size == len(df):
            df["time"] = pd.to_datetime(np.ravel(ds[tn].values), errors="coerce")
        return df
    return pd.read_csv(f, sep=None, engine="python")


def load_floats(dir_floats):
    """{float id: (lon, lat)}; one file per float, or one file with an id column."""
    floats = {}
    files = sorted({f for p in FLOAT_PATTERNS for f in Path(dir_floats).glob(p)},
                   key=lambda f: natural_key(f.name))
    for f in files:
        try:
            df = _read_any(f)
            clon, clat = lonlat_cols(df)
        except Exception as e:                              # noqa: BLE001
            print(f"  ! float file {f.name} skipped ({e})")
            continue
        tc, ic = time_col(df), id_col(df)
        if tc is not None:
            df = df.assign(_t=pd.to_datetime(df[tc], errors="coerce")).sort_values("_t")
        groups = df.groupby(ic, sort=False) if ic else [(f.stem, df)]
        for fid, g in groups:
            lo, la = clean_lonlat(g[clon], g[clat])
            if np.isfinite(lo).any():
                floats[str(fid)] = (lo[np.isfinite(lo)], la[np.isfinite(la)])
    print(f"  floats: {len(floats)}")
    return floats


def load_drifters(dir_drifters=DIR_DRIFTERS):
    """{drifter id: DataFrame(time, lon, lat)} between FIRST_LAUNCH and DRIFTER_END."""
    dfs = []
    files = sorted(Path(dir_drifters).glob("*.csv"))
    for k, f in enumerate(files):
        head = pd.read_csv(f, nrows=0)                      # read only the needed columns
        keep = [c for c in (*lonlat_cols(head), time_col(head), id_col(head)) if c]
        d = pd.read_csv(f, usecols=keep).drop_duplicates()
        d["_file"] = f.stem
        dfs.append(d)
        if (k + 1) % 200 == 0:
            tick(f"drifter files read: {k + 1}/{len(files)}")
    df = pd.concat(dfs, ignore_index=True).drop_duplicates()
    clon, clat = lonlat_cols(df)
    tc = time_col(df)
    if tc is None:
        raise KeyError(f"no time column in drifter files (columns: {list(df.columns)})")
    ic = id_col(df) or "_file"
    out = pd.DataFrame({"id": df[ic].astype(str),
                        "time": pd.to_datetime(df[tc], errors="coerce", utc=True).dt.tz_localize(None)})
    out["lon"], out["lat"] = clean_lonlat(df[clon], df[clat])
    out = out.dropna()
    out = out[out["time"] >= pd.Timestamp(FIRST_LAUNCH)]
    if DRIFTER_END:
        out = out[out["time"] <= pd.Timestamp(DRIFTER_END) + pd.Timedelta("1D")]
    out = out.drop_duplicates(["id", "time"]).sort_values(["id", "time"])
    drifters = {k: g.reset_index(drop=True) for k, g in out.groupby("id")}
    print(f"  drifters: {len(drifters)}")
    return drifters


def drifter_segments(d):
    """Line segments (n, 2, 2) of one drifter and the speed (m/s) along each segment."""
    t = d["time"].values.astype("datetime64[s]").astype(float)
    lon, lat = d["lon"].values, d["lat"].values
    if len(t) < 2:
        return np.empty((0, 2, 2)), np.empty(0)
    dt = np.diff(t)
    dx = np.deg2rad(np.diff(lon)) * 6371e3 * np.cos(np.deg2rad(0.5 * (lat[1:] + lat[:-1])))
    dy = np.deg2rad(np.diff(lat)) * 6371e3
    with np.errstate(divide="ignore", invalid="ignore"):
        spd = np.hypot(dx, dy) / dt
    ok = (dt > 0) & (dt <= MAX_DT_H * 3600) & (spd <= MAX_SPEED)
    spd = np.where(ok, spd, np.nan)
    if SPEED_SMOOTH_N > 1:
        spd = pd.Series(spd).rolling(SPEED_SMOOTH_N, center=True, min_periods=1).mean().values
    segs = np.stack([np.c_[lon[:-1], lat[:-1]], np.c_[lon[1:], lat[1:]]], axis=1)
    keep = ok & np.isfinite(spd)
    return segs[keep], spd[keep]


# =============================================================================
# Plot pieces
# =============================================================================
def plot_geoms(ax, gdf, **kw):
    for geom in gdf.geometry:
        if geom is None:
            continue
        for poly in getattr(geom, "geoms", [geom]):
            x, y = poly.exterior.xy if hasattr(poly, "exterior") else poly.xy
            ax.plot(x, y, transform=crs, **kw)


def draw_track(ax, lonlat, colour, lw=1.3, z=7):
    lon, lat = lonlat
    ax.plot(lon, lat, color=colour, lw=lw, transform=crs, zorder=z,
            solid_capstyle="round", path_effects=halo(lw + 1.6))


def draw_background(ax, extent, sst, sst_kw, eddies, grey=False):
    """SST + isobaths + TOEddies + land + grid. Returns the SST mappable."""
    ax.set_extent(extent, crs=crs)
    pc = None
    if sst is not None:
        kw = dict(sst_kw, cmap="Greys_r") if grey else sst_kw
        pc = ax.pcolormesh(sst.lon, sst.lat, sst.values, shading="auto", transform=crs,
                           zorder=1, rasterized=True, **kw)
    if SHOW["bathy"] and FNAME_BATHY.exists():
        b = xr.open_dataset(FNAME_BATHY)["topo"]
        b = b.sel(lon=slice(extent[0] - 1, extent[1] + 1), lat=slice(extent[2] - 1, extent[3] + 1))
        ax.contour(b.lon, b.lat, b.values, levels=sorted(BATHY_LEVELS), colors="0.85",
                   linewidths=0.5, alpha=0.7, transform=crs, zorder=2)
    for kind, col in (("anti", AC_COLOUR), ("cyc", C_COLOUR)):
        e = (eddies or {}).get(kind, {})
        if "maxspeed" in e:
            plot_geoms(ax, e["maxspeed"], color=col, lw=1.6, zorder=5)
        if SHOW["eddy_outer"] and "outcontour" in e:
            plot_geoms(ax, e["outcontour"], color=col, lw=0.8, ls=(0, (4, 2)), zorder=5)
    ax.add_feature(cfeature.LAND.with_scale(LAND_RES), facecolor="0.82", edgecolor="none", zorder=4)
    ax.coastlines(resolution=LAND_RES, linewidth=0.6, color="0.25", zorder=4)
    gl = ax.gridlines(crs=crs, draw_labels=True, linewidth=0.4, color="0.4", alpha=0.35,
                      linestyle="--", zorder=3)
    gl.top_labels = gl.right_labels = False
    return pc


def auto_extent(points, pad=EXTENT_PAD):
    pts = [p for p in points if p is not None and len(p[0])]
    if not pts:
        return [8, 25, -41, -34]
    lon = np.concatenate([np.ravel(p[0]) for p in pts])
    lat = np.concatenate([np.ravel(p[1]) for p in pts])
    r = lambda v, f: f(v * 2) / 2                       # round to 0.5 deg
    return [r(np.nanmin(lon) - pad, np.floor), r(np.nanmax(lon) + pad, np.ceil),
            r(np.nanmin(lat) - pad, np.floor), r(np.nanmax(lat) + pad, np.ceil)]


def header(text):
    """Invisible handle used as a bold section title inside the legend."""
    return mlines.Line2D([], [], ls="none", label=text)


def line_handle(colour, label, lw=1.6, **kw):
    return mlines.Line2D([], [], color=colour, lw=lw, label=label,
                         path_effects=kw.pop("path_effects", halo(lw + 1.6)), **kw)


def eddy_handles(eddy_key):
    if not eddy_key:
        return []
    h = [header(f"TOEddies, {DATE}"),
         mlines.Line2D([], [], color=AC_COLOUR, lw=1.6, label="Anticyclone (max speed)"),
         mlines.Line2D([], [], color=C_COLOUR, lw=1.6, label="Cyclone (max speed)",
                       path_effects=halo(4, "0.55"))]
    if SHOW["eddy_outer"]:
        h.append(mlines.Line2D([], [], color="0.3", lw=0.8, ls=(0, (4, 2)),
                               label="Outer contour"))
    return h


def make_legend(ax, handles, ncol, **kw):
    ax.axis("off")
    leg = ax.legend(handles=handles, ncol=ncol, frameon=False, fontsize=FONTSIZE - 1,
                    handlelength=2.2, columnspacing=1.6, labelspacing=0.5, **kw)
    handles_drawn = getattr(leg, "legend_handles", None) or leg.legendHandles  # mpl < 3.7
    for h, t in zip(handles_drawn, leg.get_texts()):
        if isinstance(h, mlines.Line2D) and h.get_linestyle() == "None" \
                and h.get_marker() in ("None", None, ""):
            t.set_fontweight("bold")
    return leg


# =============================================================================
# Main
# =============================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading data")

    # ---- left-panel platforms ----
    gliders = []
    if SHOW["gliders"]:
        for f in sorted(DIR_SG.glob(GLIDER_GLOB), key=lambda f: natural_key(f.name)):
            sg = re.search(r"sg\s*_?(\d+)", f.stem, re.I)
            gliders.append((f"Seaglider SG{sg.group(1)}" if sg else f.stem, read_track(f)))
        gliders.append(("SeaExplorer", read_track(SEAEXPLORER_FILE)))
        gliders = [(lab, tr) for lab, tr in gliders if tr is not None]
    wgs = []
    if SHOW["wavegliders"]:
        wgs = [(lab, read_track(f)) for f, lab in WAVEGLIDERS]
        wgs = [(lab, tr) for lab, tr in wgs if tr is not None]
    xspar = read_track(XSPAR_FILE) if SHOW["xspar"] else None
    floats = load_floats(DIR_FLOATS) if SHOW["floats"] and DIR_FLOATS.exists() else {}
    ctd = load_ctd(DIR_CTD_META) if SHOW["ctd"] and DIR_CTD_META.exists() else None
    tick("gliders, wave gliders, XSPAR, floats, CTD loaded")
    mvp = load_mvp(DIR_MVP) if SHOW["mvp"] and DIR_MVP.exists() else {}
    tick("MVP loaded")

    # ---- right-panel drifters ----
    drifters = {}
    if SHOW["drifters"]:
        try:
            drifters = load_drifters()
        except Exception as e:                              # noqa: BLE001
            print(f"  ! drifters skipped ({e})")

    tick("drifters loaded")
    pts = [tr for _, tr in gliders + wgs] + [xspar, ctd] + list(mvp.values()) \
        + list(floats.values()) + [(d["lon"].values, d["lat"].values) for d in drifters.values()]
    extent = EXTENT or auto_extent(pts)
    print(f"  map extent: {extent}")

    # ---- shared background data ----
    sst, sst_kw = None, {}
    if SHOW["sst"]:
        sst = load_sst(SST_FILE, DATE).sel(lon=slice(extent[0] - 1, extent[1] + 1),
                                          lat=slice(extent[2] - 1, extent[3] + 1))
        if SST_RANGE:
            vmin, vmax = SST_RANGE
        else:
            vmin, vmax = np.nanpercentile(sst.values, [2, 98])
            vmin, vmax = np.floor(vmin * 2) / 2, np.ceil(vmax * 2) / 2
        sst_kw = dict(cmap=SST_CMAP, vmin=vmin, vmax=vmax)
    eddies, eddy_key = None, None
    if SHOW["eddies"]:
        files, eddy_key = get_eddy_files(DATE)
        eddies = load_eddies(files) if files else None
    tick("SST + eddies loaded")

    # ---- layout: 2 maps on top, 2 legend boxes below ----
    proj = ccrs.Mercator(central_longitude=0.5 * (extent[0] + extent[1]))
    fig = plt.figure(figsize=FIGSIZE, constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.22])
    axa = fig.add_subplot(gs[0, 0], projection=proj)
    axb = fig.add_subplot(gs[0, 1], projection=proj)
    axla, axlb = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    pc = draw_background(axa, extent, sst, sst_kw, eddies)
    draw_background(axb, extent, sst, sst_kw, eddies, grey=DRIFTER_PANEL_SST == "grey")
    tick("backgrounds drawn")

    # ================= (a) platforms =================
    for (_, tr), col in zip(gliders, GLIDER_COLOURS):
        draw_track(axa, tr, col)
    for (_, tr), col in zip(wgs, WAVEGLIDER_COLOURS):
        draw_track(axa, tr, col)
    if xspar is not None:
        draw_track(axa, xspar, XSPAR_COLOUR)
    for lo, la in floats.values():
        draw_track(axa, (lo, la), FLOAT_COLOUR, lw=1.1)
        axa.scatter(lo[0], la[0], s=40, marker="^", facecolor=FLOAT_COLOUR, edgecolor="white",
                    lw=0.8, transform=crs, zorder=9)
    for name, (lo, la) in mvp.items():
        draw_track(axa, (lo, la), MVP_COLOUR, lw=2.2, z=8)
        if LABEL_MVP and np.isfinite(lo).any():
            i = np.flatnonzero(np.isfinite(lo))[0]
            num = re.search(r"\d+", name)
            axa.annotate(num.group(0) if num else name, (lo[i], la[i]),
                         xycoords=crs._as_mpl_transform(axa), xytext=(-3, 3),
                         textcoords="offset points", ha="right", va="bottom", fontsize=7,
                         fontweight="bold", zorder=12, annotation_clip=True,
                         path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
    if ctd is not None:
        axa.scatter(*ctd, s=22, marker="o", facecolor=CTD_FACE, edgecolor="k", lw=0.8,
                    transform=crs, zorder=10)
    axa.set_title("(a) Ship-based and autonomous platforms", loc="left")

    hs = [header("Ship-based")]
    if ctd is not None:
        hs.append(mlines.Line2D([], [], ls="none", marker="o", ms=5, mfc=CTD_FACE, mec="k",
                                label=f"CTD stations (n = {np.isfinite(ctd[0]).sum()})"))
    if mvp:
        hs.append(line_handle(MVP_COLOUR, f"MVP transects (n = {len(mvp)})"
                              + (", numbered" if LABEL_MVP else ""), lw=2.2))
    hs.append(header("Autonomous"))
    hs += [line_handle(c, lab) for (lab, _), c in zip(gliders, GLIDER_COLOURS)]
    hs += [line_handle(c, lab) for (lab, _), c in zip(wgs, WAVEGLIDER_COLOURS)]
    if xspar is not None:
        hs.append(line_handle(XSPAR_COLOUR, "XSPAR (Xeos)"))
    if floats:
        hs.append(line_handle(FLOAT_COLOUR, f"Floats (n = {len(floats)}), ▲ deployment",
                              lw=1.1))
    hs += eddy_handles(eddy_key)
    make_legend(axla, hs, ncol=3, loc="upper center")

    # ================= (b) drifters =================
    segs, spd = [], []
    for d in drifters.values():
        s, v = drifter_segments(d)
        segs.append(s)
        spd.append(v)
    hb = [header("Surface drifters")]
    if drifters:
        segs, spd = np.concatenate(segs), np.concatenate(spd)
        smin, smax = SPEED_LIM or (0.0, max(0.2, np.ceil(np.nanpercentile(spd, 98) * 10) / 10))
        lc = LineCollection(segs, array=spd, cmap=SPEED_CMAP, norm=Normalize(smin, smax),
                            linewidths=1.5, transform=crs, zorder=8, capstyle="round",
                            path_effects=halo(3.2, "0.1"))
        axb.add_collection(lc)
        dep = np.array([(d["lon"].iloc[0], d["lat"].iloc[0]) for d in drifters.values()])
        axb.scatter(dep[:, 0], dep[:, 1], s=80, marker="*", facecolor=DRIFTER_FACE,
                    edgecolor="k", lw=0.6, transform=crs, zorder=11)
        hb.append(mlines.Line2D([], [], ls="none", marker="*", ms=11, mfc=DRIFTER_FACE, mec="k",
                                label=f"Deployment (n = {len(drifters)})"))
        # speed colour bar inside the legend box, so both maps keep the same size
        cax = axlb.inset_axes([0.55, 0.62, 0.4, 0.14])
        cbs = fig.colorbar(lc, cax=cax, orientation="horizontal", extend="max")
        cbs.set_label("Drifter speed (m s$^{-1}$)")
    end = pd.Timestamp(DRIFTER_END).strftime("%d %b") if DRIFTER_END else "end of record"
    axb.set_title(f"(b) Surface drifters, coloured by speed "
                  f"({pd.Timestamp(FIRST_LAUNCH):%d %b} to {end})", loc="left")
    make_legend(axlb, hb, ncol=2, loc="upper left", bbox_to_anchor=(0.0, 1.0))

    # ---- shared SST colour bar under both maps ----
    if pc is not None:
        cb = fig.colorbar(pc, ax=[axa, axb], orientation="horizontal", shrink=0.45, pad=0.01,
                          aspect=50, extend="both")
        cb.set_label(f"Sea surface temperature (°C), CLS 7-day composite, {DATE}"
                     + (" (grey in b)" if DRIFTER_PANEL_SST == "grey" else ""))
    eddy_txt = {"TOEDDIES_better_ADT": "TOEddies (better ADT)", "TOEDDIES": "TOEddies"}.get(eddy_key, "eddies")
    fig.suptitle(f"WHIRLS sampling overview: SST and {eddy_txt} on {DATE}; "
                 "platform positions for the whole cruise", fontsize=FONTSIZE + 3)

    tick("everything plotted (now rendering)")
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{OUT_NAME}.{ext}", dpi=300, bbox_inches="tight")
        tick(f"saved .{ext}")
    print(f"Saved {OUT_DIR / OUT_NAME}.png / .pdf")
    plt.show()


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()