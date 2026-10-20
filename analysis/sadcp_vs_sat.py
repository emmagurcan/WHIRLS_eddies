# %% [markdown]
# # WHIRLS 2026 — SADCP vs altimetry-derived surface geostrophic velocity
#
# Along-track comparison in the spirit of Subirade et al. (2023, Remote Sensing 15, 1897),
# Figures 7-8 and Table 2:
#   1. Surface geostrophic velocity computed from daily ADT maps
#        ug = -(g/f) d(ADT)/dy ,  vg = (g/f) d(ADT)/dx
#      and interpolated (space + time) onto the ship track.
#   2. SADCP (attitude-corrected OS150 for the upper ~200 m, OS38 deeper), QC'd with
#      CAS_CURRENT_FLAG == 1, averaged in time windows (default 1 h).
#   3. Metrics: speed difference |U_sadcp| - |U_geo|, angle difference (SADCP relative to
#      geostrophy, positive = SADCP rotated anticlockwise), Pearson correlation vs depth,
#      complex (Kundu 1976) vector correlation and mean veering angle.
#   4. Figures: Fig.7-like section panel, Fig.8-like scatter, Table-2-like statistics,
#      and an animation where SADCP (red) and geostrophic (blue) arrows accumulate during a
#      day and are cleared when the next day starts.
#
# Run as a script (python sadcp_vs_geostrophy.py) or cell by cell in VS Code / Jupyter.

# %% Imports
import os
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib

if os.environ.get("MPLBACKEND") is None and not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines
import cartopy.crs as ccrs
import cartopy.feature as cfeature

try:
    import geopandas as gpd
except ImportError:  # eddy contours are optional
    gpd = None

# %% Configuration — edit paths here
BASE = Path(os.environ.get("WHIRLS_BASE", "/home/geosciences/Documents/Gurcan"))
DIR_SADCP = BASE / "data/SADCP"
DIR_SSH = BASE / "data/SSH"
DIR_EDDIES = BASE / "data/TOEDDIES/TOEDDIES_better_ADT"
FNAME_BATHY = BASE / "WHIRLS_eddies/REP/data/bathymetry/etopo2.nc"
DIR_OUT = BASE / "WHIRLS_eddies/figures/sadcp_vs_geostrophy"

# README: the *_TOT_osite files hold the most recent (full-cruise) data; the MD2_* files are
# attitude-corrected but may cover only part of the cruise.
# SADCP_SOURCE: "TOT"      -> TOT files only (full cruise)
#               "MD2"      -> attitude-corrected files only
#               "MD2+TOT"  -> MD2 where available, gaps filled with TOT  (default)
SADCP_SOURCE = "MD2+TOT"
SADCP_FILES = {
    "shallow": {"MD2": DIR_SADCP / "MD2_OS150_WHIRLS-2026_120s_oceansites.nc",   # ~200 m, 8 m bins
                "TOT": DIR_SADCP / "WHIRLS_NB150_TOT_osite.nc"},
    "deep":    {"MD2": DIR_SADCP / "MD2_OS38_WHIRLS-2026_120s_oceansites.nc",    # ~1000 m, 24 m bins
                "TOT": DIR_SADCP / "WHIRLS_BB38_TOT_osite.nc"},
}

CRUISE_START = "2026-06-28"
CRUISE_END = "2026-07-27T23:59:59"

AVG_WINDOW = "1h"          # SADCP time averaging (reduces ping noise, not tides/inertial)
MIN_COVERAGE = 0.25        # keep a window if >= 25% of its expected ensembles are valid
                           # (expected count comes from each file's own sampling interval)
MIN_LAYER_FRAC = 0.5       # 0-100 m mean needs >= 50% of the bins that normally hold data
MAX_ADT_GAP = pd.Timedelta(hours=36)  # no geostrophy if the nearest ADT map is further away
Z_SURF = 30.0              # m — "near-surface" level (green arrows in the paper)
Z_LAYER = (50.0, 100.0)     # m — layer average (red arrows in the paper)
Z_SPLIT = 150.0            # m — use OS150 above, OS38 below for depth profiles
DEPTHS_SECTION = [30, 50, 75, 100, 150, 200, 300, 500]  # levels shown in panels b, d, e
MIN_SPEED_ANGLE = 0.10     # m/s — angle differences below this speed are unreliable
ARROW_EVERY = "3h"         # arrow thinning on static maps
ARROW_SCALE = 2.0          # m/s per inch of arrow
ANIM_STEP = "2h"           # one animation frame every ANIM_STEP
ANIM_LEVEL = "layer"       # "layer" (0-100 m mean) or "surf" (30 m) for the ship arrows
ANIM_FPS = 6
ADD_LAND = os.environ.get("WHIRLS_NO_LAND") is None

# Optional named sections (e.g. eddy crossings) for extra Fig.7-like panels:
SECTIONS = {
    # "eddy_A_crossing": ("2026-07-10T00:00", "2026-07-12T18:00"),
}

G = 9.81
OMEGA = 7.2921e-5
R_EARTH = 6371e3

DIR_OUT.mkdir(parents=True, exist_ok=True)


# %% Helpers: variable lookup
def _find(ds, names, required=True):
    for n in names:
        if n in ds.variables:
            return n
    low = {k.lower(): k for k in ds.variables}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    if required:
        raise KeyError(f"None of {names} found. Variables: {list(ds.variables)}")
    return None


# %% Load SADCP
def load_sadcp(path, t0=CRUISE_START, t1=CRUISE_END, window=AVG_WINDOW, verbose=True):
    """Return Dataset(time, z) with u, v, n (valid ensembles) and lat/lon(time)."""
    ds = xr.open_dataset(path)
    if verbose:
        print(f"\n--- {Path(path).name} ---")
        print({k: (ds[k].dims, ds[k].shape) for k in ds.data_vars})

    tn = _find(ds, ["JULD_ADCP", "TIME", "time", "JULD"])
    zn = _find(ds, ["DEPH", "DEPTH", "depth", "bin_depth", "PRES"])
    yn = _find(ds, ["LATITUDE", "latitude", "lat"])
    xn = _find(ds, ["LONGITUDE", "longitude", "lon"])
    un = _find(ds, ["UVEL_ADCP", "UCUR", "u"])
    vn = _find(ds, ["VVEL_ADCP", "VCUR", "v"])
    fn = _find(ds, ["CAS_CURRENT_FLAG", "UVEL_ADCP_QC", "UCUR_QC"], required=False)

    tdim = ds[tn].dims[0]
    u = ds[un].squeeze(drop=True)
    zdim = [d for d in u.dims if d != tdim][0]
    u = u.transpose(tdim, zdim).values.astype(float)
    v = ds[vn].squeeze(drop=True).transpose(tdim, zdim).values.astype(float)
    traw = ds[tn].values
    if np.issubdtype(traw.dtype, np.datetime64):
        time = pd.to_datetime(traw)
    else:  # JULD_ADCP: decimal days since 1950-01-01 (not decoded by xarray)
        time = pd.Timestamp("1950-01-01") + pd.to_timedelta(traw.astype(float), unit="D")

    if fn is not None:
        fl = ds[fn].squeeze(drop=True)
        fl = fl.transpose(tdim, zdim).values if zdim in fl.dims else fl.values[:, None]
        if np.issubdtype(fl.dtype, np.number):
            good = fl == 1
        else:  # char/bytes flags, e.g. b'1'
            good = np.char.strip(np.char.decode(fl.astype("S"), "ascii", "ignore")) == "1"
        good = np.broadcast_to(good, u.shape)
        u = np.where(good, u, np.nan)
        v = np.where(good, v, np.nan)
    else:
        print("  ! no QC flag variable found, using all data")
    bad = (np.abs(u) > 5) | (np.abs(v) > 5)
    u[bad] = np.nan
    v[bad] = np.nan

    # depth: 1D (bins) or 2D (time, bins)
    zda = ds[zn].squeeze(drop=True)
    if zda.ndim == 2:
        z2 = zda.transpose(tdim, zdim).values.astype(float)
        zgrid = np.round(np.nanmedian(z2, axis=0), 1)
        u, v = _to_grid(u, z2, zgrid), _to_grid(v, z2, zgrid)
    else:
        zgrid = zda.values.astype(float)
    zgrid = np.abs(zgrid)  # depth positive down
    zok = np.isfinite(zgrid)              # drop bins with no defined depth
    u, v, zgrid = u[:, zok], v[:, zok], zgrid[zok]

    lat = ds[yn].squeeze(drop=True).values.astype(float)
    lon = ds[xn].squeeze(drop=True).values.astype(float)
    lat = np.where(np.abs(lat) <= 90, lat, np.nan)          # undeclared fill values
    lon = np.where(np.abs(lon) <= 360, lon, np.nan)
    lon = np.where(lon > 180, lon - 360, lon)

    out = xr.Dataset(
        {"u": (("time", "z"), u), "v": (("time", "z"), v),
         "lat": ("time", lat), "lon": ("time", lon)},
        coords={"time": time, "z": zgrid},
    )
    out = out.isel(time=np.asarray(pd.notna(time)))
    out = out.sortby("time").sel(time=slice(t0, t1))
    out = out.isel(time=~pd.Index(out.time.values).duplicated())

    # expected ensembles per window from the file's own sampling interval
    dt = pd.Series(out.time.values).diff().median()
    expected = max(1.0, pd.Timedelta(window) / dt) if pd.notna(dt) and dt > pd.Timedelta(0) else 1.0
    n_min = max(1, int(np.ceil(MIN_COVERAGE * expected)))

    ok = np.isfinite(out.u)
    n = ok.resample(time=window).sum()
    avg = out.resample(time=window).mean()
    avg["n"] = n
    avg["u"] = avg.u.where(n >= n_min)
    avg["v"] = avg.v.where(n >= n_min)
    avg = avg.sortby("z")
    if verbose:
        raw_good = float(ok.mean()) * 100
        print(f"  time {pd.Timestamp(out.time.values[0]):%Y-%m-%d %H:%M} -> "
              f"{pd.Timestamp(out.time.values[-1]):%Y-%m-%d %H:%M}, sampling {dt}, "
              f"{len(out.time)} ensembles ({raw_good:.0f}% of cells pass QC)")
        print(f"  -> {len(avg.time)} windows of {window} (need >= {n_min} of ~{expected:.0f}); "
              f"z = {avg.z.values.min():.0f}-{avg.z.values.max():.0f} m ({len(avg.z)} bins)")
    ds.close()
    return avg


def load_sadcp_source(kind, source=SADCP_SOURCE):
    """Load the shallow/deep SADCP from TOT, MD2, or MD2 filled with TOT."""
    files = SADCP_FILES[kind]
    if source in ("TOT", "MD2"):
        return load_sadcp(files[source]) if files[source].exists() else None
    md2 = load_sadcp(files["MD2"]) if files["MD2"].exists() else None
    tot = load_sadcp(files["TOT"]) if files["TOT"].exists() else None
    if md2 is None or tot is None:
        return md2 if tot is None else tot
    md2 = md2.interp(z=tot.z)                 # put MD2 on the TOT bin grid
    merged = md2.combine_first(tot)           # MD2 wins, TOT fills the gaps
    used_md2 = np.isfinite(md2.u.reindex(time=merged.time)).any("z")
    print(f"  {kind}: MD2 used for {int(used_md2.sum())} windows, "
          f"TOT fills {int((np.isfinite(merged.u).any('z') & ~used_md2).sum())}")
    return merged


def _to_grid(a, z2, zgrid):
    out = np.full((a.shape[0], len(zgrid)), np.nan)
    for i in range(a.shape[0]):
        m = np.isfinite(a[i]) & np.isfinite(z2[i])
        if m.sum() >= 2:
            zi, ai = z2[i, m], a[i, m]
            o = np.argsort(zi)
            out[i] = np.interp(zgrid, zi[o], ai[o], left=np.nan, right=np.nan)
    return out


def level(sad, z0):
    """SADCP u, v at the bin closest to z0."""
    s = sad.sel(z=z0, method="nearest").drop_vars("z")
    return s.u, s.v


def layer_mean(sad, z1, z2, min_frac=MIN_LAYER_FRAC):
    """Depth mean over [z1, z2]. Bins that are almost never valid (e.g. the blanked top
    bins of a ship ADCP) are ignored so they don't make every profile fail."""
    s = sad.sel(z=slice(z1, z2))
    usable = np.isfinite(s.u).mean("time") >= 0.10
    s = s.sel(z=usable)
    if s.sizes["z"] == 0:
        nan = xr.full_like(sad.lat, np.nan)
        return nan, nan
    frac = np.isfinite(s.u).sum("z") / s.sizes["z"]
    return s.u.mean("z").where(frac >= min_frac), s.v.mean("z").where(frac >= min_frac)


# %% Load ADT and compute geostrophy
def load_adt(dir_ssh=DIR_SSH, t0=CRUISE_START, t1=CRUISE_END, bbox=None, pad_days=1):
    files = sorted(glob.glob(str(Path(dir_ssh) / "*.nc")))
    days = pd.date_range(pd.Timestamp(t0) - pd.Timedelta(days=pad_days),
                         pd.Timestamp(t1) + pd.Timedelta(days=pad_days), freq="D")
    tags = {d.strftime("%Y-%m-%d") for d in days} | {d.strftime("%Y%m%d") for d in days}
    picked = [f for f in files if any(t in Path(f).name for t in tags)] or files
    print(f"ADT: {len(picked)} file(s) from {dir_ssh}")

    dsets = []
    for f in picked:
        d = xr.open_dataset(f)
        ren = {k: v for k, v in {"latitude": "lat", "longitude": "lon"}.items() if k in d.dims}
        d = d.rename(ren)
        an = _find(d, ["adt", "ADT", "zos", "ssh"])
        keep = [an] + [k for k in ("ugos", "vgos") if k in d]
        d = d[keep].rename({an: "adt"})
        # daily files (e.g. 2026-07-08_cmems_data.nc): take the date from the file name,
        # as the old notebook did by dropping the file's own time variable
        if d.sizes.get("time", 1) == 1:
            if "time" in d.dims:
                d = d.isel(time=0, drop=True)
            elif "time" in d.coords:
                d = d.drop_vars("time")
            d = d.expand_dims(time=[_date_from_name(f)])
        dsets.append(d.load())
        d.close()
    ds = xr.concat(dsets, "time")
    ds = ds.assign_coords(time=np.asarray(ds.time.values, dtype="datetime64[ns]")).sortby("time")
    ds = ds.isel(time=~pd.Index(ds.time.values).duplicated())
    ds = ds.assign_coords(lon=xr.where(ds.lon > 180, ds.lon - 360, ds.lon)).sortby(["lat", "lon"])
    if bbox is not None:
        x0, x1, y0, y1 = bbox
        ds = ds.sel(lon=slice(x0, x1), lat=slice(y0, y1))
    return ds.sel(time=slice(days[0], days[-1] + pd.Timedelta(hours=23)))


def _date_from_name(f):
    m = re.search(r"(\d{4})-?(\d{2})-?(\d{2})", Path(f).name)
    return pd.Timestamp(f"{m.group(1)}-{m.group(2)}-{m.group(3)}")


def geostrophy(adt):
    """Surface geostrophic velocity (m/s) from ADT (m) on a regular lon/lat grid."""
    lat = np.deg2rad(adt.lat)
    f = 2 * OMEGA * np.sin(lat)
    deg = np.pi / 180
    dadt_dy = adt.differentiate("lat") / (R_EARTH * deg)
    dadt_dx = adt.differentiate("lon") / (R_EARTH * np.cos(lat) * deg)
    ug = -(G / f) * dadt_dy
    vg = (G / f) * dadt_dx
    return ug.rename("ug"), vg.rename("vg")


def interp_to_track(geo, time, lat, lon):
    """Linear interpolation in time, lat and lon onto ship positions.
    Points whose nearest ADT map is more than MAX_ADT_GAP away are set to NaN (no silent
    extrapolation from an old map)."""
    # Time is interpolated as plain float days: xarray's own datetime interpolation can
    # mix time units (s/us vs ns) between the maps and the query and return all-NaN.
    epoch = np.datetime64("2000-01-01T00:00:00", "ns")
    day = np.timedelta64(86400, "s")
    t = (np.asarray(time, dtype="datetime64[ns]") - epoch) / day
    tg = (np.asarray(geo.time.values, dtype="datetime64[ns]") - epoch) / day
    geo_num = geo.assign_coords(time=tg)
    tq = np.clip(t, tg.min(), tg.max())
    out = geo_num.interp(time=xr.DataArray(tq, dims="obs"),
                         lat=xr.DataArray(np.asarray(lat, dtype=float), dims="obs"),
                         lon=xr.DataArray(np.asarray(lon, dtype=float), dims="obs"))
    idx = np.clip(np.searchsorted(tg, t), 1, len(tg) - 1) if len(tg) > 1 else np.zeros(len(t), int)
    gap = np.minimum(np.abs(t - tg[idx]), np.abs(t - tg[np.maximum(idx - 1, 0)]))
    far = xr.DataArray(gap > MAX_ADT_GAP / pd.Timedelta(days=1), dims="obs")
    out = out.where(~far)
    if not np.isfinite(out["ug"].values).any():
        print("  ! geostrophy is NaN on the whole track. Check:")
        print(f"    track lon {np.nanmin(lon):.2f}..{np.nanmax(lon):.2f}, "
              f"lat {np.nanmin(lat):.2f}..{np.nanmax(lat):.2f}")
        print(f"    grid  lon {float(geo.lon.min()):.2f}..{float(geo.lon.max()):.2f}, "
              f"lat {float(geo.lat.min()):.2f}..{float(geo.lat.max()):.2f}, "
              f"NaN fraction of ug grid {float(np.isnan(geo.ug).mean()):.2f}")
    return out


def adt_for_day(geo, day):
    """ADT map for a day, or None if no map within MAX_ADT_GAP."""
    a = geo.adt.sel(time=day, method="nearest")
    if abs(pd.Timestamp(a.time.values) - pd.Timestamp(day)) > MAX_ADT_GAP:
        return None
    return a


def coverage_report(sh, geo, fname="coverage_by_day.csv"):
    """Per-day count of hourly windows that are valid at each step, to locate NaNs."""
    df = pd.DataFrame({
        "position": np.isfinite(sh.lat.values) & np.isfinite(sh.lon.values),
        "sadcp_any_bin": np.isfinite(sh.u.values).any(axis=1),
        "sadcp_30m": np.isfinite(sh.u30.values),
        "sadcp_0_100m": np.isfinite(sh.u100.values),
        "geostrophy": np.isfinite(sh.ug.values),
    }, index=pd.to_datetime(sh.time.values))
    df["both_0_100m"] = df.sadcp_0_100m & df.geostrophy
    rep = df.groupby(df.index.floor("D")).sum().astype(int)
    rep.insert(0, "windows", df.groupby(df.index.floor("D")).size())
    days = pd.DatetimeIndex(pd.to_datetime(geo.time.values)).floor("D")
    rep["adt_map"] = rep.index.isin(days)
    rep.to_csv(DIR_OUT / fname)
    print("\nCoverage per day (number of hourly windows that are valid):")
    print(rep.to_string())
    return rep


# %% Metrics
def speed_angle_diff(uo, vo, ug, vg):
    """Return numpy arrays: |U_obs|-|U_geo| (m/s) and angle(U_obs) - angle(U_geo) in [-180, 180]."""
    uo, vo, ug, vg = (np.asarray(a, dtype=float) for a in (uo, vo, ug, vg))
    wo, wg = uo + 1j * vo, ug + 1j * vg
    dnorm = np.abs(wo) - np.abs(wg)                       # >0: SADCP faster
    dang = np.degrees(np.angle(wo * np.conj(wg)))         # SADCP direction - geo direction
    return dnorm, dang


def pearson(a, b, nmin=10):
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    m = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[m], b[m])[0, 1] if m.sum() >= nmin else np.nan


def kundu(uo, vo, ug, vg, nmin=10):
    """Complex vector correlation magnitude and mean veering angle (SADCP rel. to geo)."""
    w1, w2 = np.asarray(ug + 1j * vg).ravel(), np.asarray(uo + 1j * vo).ravel()
    m = np.isfinite(w1) & np.isfinite(w2)
    if m.sum() < nmin:
        return np.nan, np.nan
    w1, w2 = w1[m] - w1[m].mean(), w2[m] - w2[m].mean()
    c = np.mean(np.conj(w1) * w2) / np.sqrt(np.mean(np.abs(w1) ** 2) * np.mean(np.abs(w2) ** 2))
    return np.abs(c), np.degrees(np.angle(c))


def along_track_km(lat, lon):
    lat, lon = np.deg2rad(np.asarray(lat)), np.deg2rad(np.asarray(lon))
    ok = np.isfinite(lat) & np.isfinite(lon)
    d = np.zeros(len(lat))
    idx = np.where(ok)[0]
    if len(idx) > 1:
        la, lo = lat[idx], lon[idx]
        a = (np.sin(np.diff(la) / 2) ** 2
             + np.cos(la[:-1]) * np.cos(la[1:]) * np.sin(np.diff(lo) / 2) ** 2)
        step = 2 * R_EARTH * np.arcsin(np.sqrt(a)) / 1e3
        d[idx] = np.concatenate([[0], np.cumsum(step)])
        d[~ok] = np.nan
    return d


# %% TOEddies contours
_EDDY_CACHE = {}


def load_eddies_for_day(day):
    """{('anticyclones'|'cyclones', 'maxspeed'|'outcontour'|'centers'): GeoDataFrame|None}"""
    tag = pd.Timestamp(day).strftime("%Y%m%d")
    if tag in _EDDY_CACHE:
        return _EDDY_CACHE[tag]
    out = {}
    for kind in ("anticyclones", "cyclones"):
        for part in ("maxspeed", "outcontour", "centers"):
            f = DIR_EDDIES / f"{tag}_{kind}_{part}.geojson"
            out[(kind, part)] = gpd.read_file(f) if (gpd is not None and f.exists()) else None
    _EDDY_CACHE[tag] = out
    return out


def _geom_xy(geom):
    gt = geom.geom_type
    if gt == "Polygon":
        return [np.asarray(geom.exterior.coords)]
    if gt in ("LineString", "LinearRing"):
        return [np.asarray(geom.coords)]
    if gt.startswith("Multi") or gt == "GeometryCollection":
        return [xy for g in geom.geoms for xy in _geom_xy(g)]
    return []


_SEG_CACHE = {}


def eddy_segments(day, extent, pad=1.0):
    """Contour vertices of the eddies that touch the map, as {(kind, part): [Nx2 arrays]}.
    Filtering to the map is what keeps the animation fast: the TOEddies files cover a much
    larger domain than the cruise area."""
    key = (pd.Timestamp(day).strftime("%Y%m%d"), tuple(np.round(extent, 2)))
    if key in _SEG_CACHE:
        return _SEG_CACHE[key]
    x0, x1, y0, y1 = extent
    out = {}
    for (kind, part), gdf in load_eddies_for_day(day).items():
        if gdf is None or part == "centers" or len(gdf) == 0:
            continue
        b = gdf.geometry.bounds
        keep = ((b.maxx >= x0 - pad) & (b.minx <= x1 + pad)
                & (b.maxy >= y0 - pad) & (b.miny <= y1 + pad)).values
        segs = []
        for geom in gdf.geometry.values[keep]:
            if geom is None or geom.is_empty:
                continue
            for xy in _geom_xy(geom):
                xy = xy[:, :2][np.isfinite(xy[:, :2]).all(axis=1)]   # drop NaN vertices
                if len(xy) > 1:
                    segs.append(xy)
        out[(kind, part)] = segs
    _SEG_CACHE[key] = out
    return out


def draw_eddies(ax, day, extent):
    """Max-speed (solid, red AE / blue CE) and outermost (dotted) contours, one collection
    per type. Returns the artists so they can be removed."""
    from matplotlib.collections import LineCollection
    colors = {"anticyclones": "tab:red", "cyclones": "tab:blue"}
    arts = []
    for (kind, part), segs in eddy_segments(day, extent).items():
        if not segs:
            continue
        if part == "maxspeed":
            lc = LineCollection(segs, colors=colors[kind], linewidths=1.4)
        else:
            lc = LineCollection(segs, colors="k", linewidths=0.8, linestyles=":")
        lc.set_transform(ccrs.PlateCarree()._as_mpl_transform(ax))
        lc.set_zorder(3)
        ax.add_collection(lc, autolim=False)
        arts.append(lc)
    return arts


# %% Base map
def base_map(ax, extent):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    if ADD_LAND:
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "10m",
                                                    edgecolor="k", facecolor="0.6"), zorder=4)
    if FNAME_BATHY.exists():
        b = xr.open_dataset(FNAME_BATHY)
        x0, x1, y0, y1 = extent
        topo = b["topo"].sel(lon=slice(x0, x1), lat=slice(y0, y1))
        ax.contour(topo.lon, topo.lat, topo, levels=[-4000, -3000, -2000, -1000],
                   colors="0.55", linewidths=0.4, transform=ccrs.PlateCarree(), zorder=1)
    gl = ax.gridlines(draw_labels=True, lw=0.3, color="0.6", alpha=0.6)
    gl.top_labels = gl.right_labels = False


def remove_contour(cs):
    try:
        cs.remove()
    except Exception:
        for c in getattr(cs, "collections", []):
            c.remove()


# %% ============================ MAIN PIPELINE ============================
if __name__ == "__main__" or "get_ipython" in globals():
    # ---- SADCP
    sh = load_sadcp_source("shallow")
    dp = load_sadcp_source("deep")

    # ---- ADT over the cruise area
    x0 = float(np.nanmin(sh.lon)) - 1.5
    x1 = float(np.nanmax(sh.lon)) + 1.5
    y0 = float(np.nanmin(sh.lat)) - 1.5
    y1 = float(np.nanmax(sh.lat)) + 1.5
    EXTENT = [x0 + 0.5, x1 - 0.5, y0 + 0.5, y1 - 0.5]
    adt_ds = load_adt(bbox=(x0, x1, y0, y1))
    adt_days = pd.DatetimeIndex(pd.to_datetime(adt_ds.time.values)).floor("D")
    missing = pd.date_range(CRUISE_START, pd.Timestamp(CRUISE_END).floor("D")).difference(adt_days)
    print(f"ADT maps: {adt_days.min():%Y-%m-%d} -> {adt_days.max():%Y-%m-%d} ({len(adt_days)} days)")
    if len(missing):
        print(f"  ! no ADT file for: {', '.join(d.strftime('%m-%d') for d in missing)}")
    ug_grid, vg_grid = geostrophy(adt_ds.adt)
    geo = xr.Dataset({"ug": ug_grid, "vg": vg_grid, "adt": adt_ds.adt})
    if "ugos" in adt_ds:  # sanity check against the product's own geostrophy
        rms = float(np.sqrt(((ug_grid - adt_ds.ugos) ** 2 + (vg_grid - adt_ds.vgos) ** 2).mean()))
        print(f"RMS difference own geostrophy vs product ugos/vgos: {rms:.3f} m/s")

    # ---- Geostrophy on the track (per SADCP window)
    def match(sad):
        g = interp_to_track(geo, sad.time.values, sad.lat.values, sad.lon.values)
        sad = sad.copy()
        for k in ("ug", "vg", "adt"):
            sad[k] = ("time", g[k].values)
        sad["dist"] = ("time", along_track_km(sad.lat.values, sad.lon.values))
        return sad

    sh = match(sh)
    if dp is not None:
        dp = match(dp)

    # ---- Reference SADCP velocities (paper: 30 m and 0-100 m mean)
    sh["u30"], sh["v30"] = level(sh, Z_SURF)
    sh["u100"], sh["v100"] = layer_mean(sh, *Z_LAYER)
    for tag in ("30", "100"):
        dn_, da_ = speed_angle_diff(sh[f"u{tag}"], sh[f"v{tag}"], sh.ug, sh.vg)
        sh[f"dnorm{tag}"] = ("time", dn_)
        sh[f"dang{tag}"] = ("time", da_)
        sh[f"spd{tag}"] = ("time", np.hypot(sh[f"u{tag}"].values, sh[f"v{tag}"].values))
    sh["spdg"] = ("time", np.hypot(sh.ug.values, sh.vg.values))
    coverage = coverage_report(sh, geo)

    # ---- Export the matched table
    cols = ["lat", "lon", "dist", "adt", "ug", "vg", "spdg",
            "u30", "v30", "spd30", "dnorm30", "dang30",
            "u100", "v100", "spd100", "dnorm100", "dang100"]
    table = sh[cols].to_dataframe()
    table = table.dropna(subset=["u100", "u30"], how="all")
    table.to_csv(DIR_OUT / "sadcp_vs_geostrophy_track.csv")
    print(f"Matched track table: {len(table)} rows -> {DIR_OUT/'sadcp_vs_geostrophy_track.csv'}")


# %% Depth profiles of correlation / differences
def profile_stats(sad, zmax=None):
    rows = []
    for z in sad.z.values:
        if zmax is not None and z > zmax:
            continue
        uo, vo = sad.u.sel(z=z).values, sad.v.sel(z=z).values
        ug, vg = sad.ug.values, sad.vg.values
        dn, da = speed_angle_diff(uo, vo, ug, vg)
        fast = np.hypot(uo, vo) >= MIN_SPEED_ANGLE
        rho, th = kundu(uo, vo, ug, vg)
        rows.append(dict(z=z, n=int(np.isfinite(dn).sum()),
                         r_u=pearson(uo, ug), r_v=pearson(vo, vg),
                         r_speed=pearson(np.hypot(uo, vo), np.hypot(ug, vg)),
                         kundu_rho=rho, kundu_theta=th,
                         dnorm_mean=np.nanmean(dn), dnorm_std=np.nanstd(dn),
                         dang_mean=np.nanmean(da[fast]) if fast.any() else np.nan,
                         dang_std=np.nanstd(da[fast]) if fast.any() else np.nan))
    return pd.DataFrame(rows)


def get_level_series(z, sh, dp):
    src = sh if (z <= Z_SPLIT or dp is None) else dp
    s = src.sel(z=z, method="nearest")
    return src, s.u, s.v, float(s.z)


# %% Fig. 7-like panel
def plot_fig7(sh, dp, geo, t0, t1, fname, title="", extent=None):
    s = sh.sel(time=slice(t0, t1))
    d = dp.sel(time=slice(t0, t1)) if dp is not None else None
    if s.sizes["time"] < 3:
        print(f"  skip {fname}: not enough data")
        return
    dist0 = np.nanmin(s.dist)
    tmid = pd.Timestamp(s.time.values[len(s.time) // 2]).floor("D")

    fig = plt.figure(figsize=(15, 11))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[1.05, 1], hspace=0.32, wspace=0.18)
    ax_m = fig.add_subplot(gs[0:2, 0], projection=ccrs.PlateCarree())
    ax_c = fig.add_subplot(gs[2, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[1, 1], sharex=ax_b)
    ax_e = fig.add_subplot(gs[2, 1], sharex=ax_b)

    # (a) map
    if extent is None:
        extent = [float(np.nanmin(s.lon)) - 0.7, float(np.nanmax(s.lon)) + 0.7,
                  float(np.nanmin(s.lat)) - 0.7, float(np.nanmax(s.lat)) + 0.7]
    base_map(ax_m, extent)
    a = adt_for_day(geo, tmid)
    if a is None:
        a = geo.adt.sel(time=tmid, method="nearest")
    pm = ax_m.pcolormesh(a.lon, a.lat, a, cmap="RdYlBu_r", shading="auto",
                         transform=ccrs.PlateCarree(), zorder=0)
    plt.colorbar(pm, ax=ax_m, shrink=0.7, pad=0.02, label=f"ADT (m), {tmid:%Y-%m-%d}")
    draw_eddies(ax_m, tmid, extent)
    ax_m.plot(s.lon, s.lat, "k-", lw=0.8, transform=ccrs.PlateCarree(), zorder=5)
    th = s.resample(time=ARROW_EVERY).first()
    kw = dict(transform=ccrs.PlateCarree(), scale=ARROW_SCALE, scale_units="inches",
              width=0.0028, zorder=6)
    X, Y = th.lon.values, th.lat.values
    ax_m.quiver(X, Y, th.ug.values, th.vg.values, color="tab:blue", **kw)
    ax_m.quiver(X, Y, th.u100.values, th.v100.values, color="tab:red", **kw)
    q = ax_m.quiver(X, Y, th.u30.values, th.v30.values, color="tab:green", **kw)
    ax_m.quiverkey(q, 0.86, 1.03, 0.5, "0.5 m/s", labelpos="E", color="k")
    ax_m.legend(handles=[mlines.Line2D([], [], color=c, lw=2, label=l) for c, l in
                         [("tab:blue", "Geostrophic (ADT)"),
                          ("tab:red", f"SADCP {Z_LAYER[0]:.0f}-{Z_LAYER[1]:.0f} m mean"),
                          ("tab:green", f"SADCP {Z_SURF:.0f} m")]],
                loc="upper left", fontsize=8, framealpha=0.9).set_zorder(10)
    ax_m.set_title("(a) " + title, loc="left")

    x = s.dist.values - dist0
    cmap = plt.get_cmap("jet")
    cols = [cmap(i / max(1, len(DEPTHS_SECTION) - 1)) for i in range(len(DEPTHS_SECTION))]

    # (b) speeds, (d) norm diff, (e) angle diff
    ax_b.plot(x, s.spdg, "k:", lw=2, label="Geostrophic")
    for z, c in zip(DEPTHS_SECTION, cols):
        src, uo, vo, zz = get_level_series(z, s, d)
        if zz > z + 30 or zz < z - 30:
            continue
        xs = src.dist.values - dist0
        sp = np.hypot(uo.values, vo.values)
        dn, da = speed_angle_diff(uo, vo, src.ug, src.vg)
        da = np.where(sp >= MIN_SPEED_ANGLE, da, np.nan)
        ax_b.plot(xs, sp, color=c, lw=1, label=f"{zz:.0f} m")
        ax_d.plot(xs, dn, color=c, lw=1)
        ax_e.plot(xs, da, color=c, lw=0.8, marker=".", ms=2)
    ax_b.set_ylabel("Speed (m/s)")
    ax_b.legend(ncol=3, fontsize=7)
    ax_b.set_title("(b) Speed along track", loc="left")
    ax_d.axhline(0, color="k", lw=0.6)
    ax_d.set_ylabel("|U$_{SADCP}$| − |U$_{geo}$| (m/s)")
    ax_d.set_title("(d) Speed difference", loc="left")
    ax_e.axhline(0, color="k", lw=0.6)
    ax_e.set_ylim(-180, 180)
    ax_e.set_yticks(np.arange(-180, 181, 90))
    ax_e.set_ylabel("Angle SADCP − geo (°)")
    ax_e.set_xlabel("Along-track distance (km)")
    ax_e.set_title(f"(e) Angle difference (|U$_{{SADCP}}$| ≥ {MIN_SPEED_ANGLE} m/s)", loc="left")
    for axx in (ax_b, ax_d):
        plt.setp(axx.get_xticklabels(), visible=False)

    # (c) correlation vs depth
    for src, lab, ls in ((s, "OS150", "-"), (d, "OS38", "--")):
        if src is None:
            continue
        p = profile_stats(src, zmax=1000)
        ax_c.plot(p.r_u, p.z, ls, color="tab:purple", label=f"u ({lab})")
        ax_c.plot(p.r_v, p.z, ls, color="tab:orange", label=f"v ({lab})")
        ax_c.plot(p.kundu_rho, p.z, ls, color="k", label=f"vector ρ ({lab})")
    ax_c.invert_yaxis()
    ax_c.set_xlim(-0.2, 1)
    ax_c.set_xlabel("Correlation SADCP vs geostrophic")
    ax_c.set_ylabel("Depth (m)")
    ax_c.grid(alpha=0.3)
    ax_c.legend(fontsize=7, loc="upper left")
    ax_c.set_title("(c) Pearson / vector correlation", loc="left")

    fig.suptitle(f"WHIRLS SADCP vs altimetry geostrophy — "
                 f"{pd.Timestamp(t0):%d %b %H:%M} to {pd.Timestamp(t1):%d %b %Y %H:%M}", y=0.995)
    fig.savefig(DIR_OUT / fname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {DIR_OUT/fname}")


# %% Fig. 8-like scatter + Table-2-like statistics
def plot_fig8_and_stats(sh, fname="fig8_like_scatter.png"):
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    rows = []
    bins = np.arange(0, 1.6, 0.05)
    for tag, c, lab in (("100", "tab:red", f"{Z_LAYER[0]:.0f}-{Z_LAYER[1]:.0f} m mean"),
                        ("30", "tab:green", f"{Z_SURF:.0f} m")):
        sp = sh[f"spd{tag}"].values
        dn = sh[f"dnorm{tag}"].values
        da = sh[f"dang{tag}"].values
        fast = sp >= MIN_SPEED_ANGLE
        axs[0].scatter(sp, dn, s=4, color=c, alpha=0.35, label=lab)
        axs[1].scatter(sp, da, s=4, color=c, alpha=0.35, label=lab)
        idx = np.digitize(sp, bins)
        bc, mn, sd, ma = [], [], [], []
        for i in range(1, len(bins)):
            m = (idx == i) & np.isfinite(dn)
            if m.sum() >= 5:
                bc.append(bins[i - 1] + 0.025)
                mn.append(dn[m].mean())
                sd.append(dn[m].std())
                ma.append(np.nanmean(da[m]))
        axs[0].errorbar(bc, mn, yerr=sd, color=c, mec="k", fmt="o-", ms=4, lw=1.2)
        axs[1].plot(bc, ma, "o-", color=c, mec="k", ms=4, lw=1.2)
        rho, theta = kundu(sh[f"u{tag}"].values, sh[f"v{tag}"].values, sh.ug.values, sh.vg.values)
        rows.append({
            "SADCP": lab, "N": int(np.isfinite(dn).sum()),
            "mean |U_sadcp| (m/s)": np.nanmean(sp),
            "mean |U_geo| (m/s)": np.nanmean(sh.spdg.where(np.isfinite(dn))),
            "speed diff mean (m/s)": np.nanmean(dn), "speed diff std (m/s)": np.nanstd(dn),
            "speed RMSD (m/s)": np.sqrt(np.nanmean(dn ** 2)),
            "vector RMSD (m/s)": float(np.sqrt(np.nanmean(
                (sh[f"u{tag}"] - sh.ug) ** 2 + (sh[f"v{tag}"] - sh.vg) ** 2))),
            f"angle diff mean (°, |U|≥{MIN_SPEED_ANGLE})": np.nanmean(da[fast]),
            f"angle diff std (°, |U|≥{MIN_SPEED_ANGLE})": np.nanstd(da[fast]),
            "median |angle diff| (°)": np.nanmedian(np.abs(da[fast])),
            "r(u)": pearson(sh[f"u{tag}"], sh.ug), "r(v)": pearson(sh[f"v{tag}"], sh.vg),
            "r(speed)": pearson(sp, sh.spdg),
            "vector corr ρ (Kundu)": rho, "veering θ (°)": theta,
            "ratio |U_geo|/|U_sadcp| (median)": np.nanmedian(
                (sh.spdg / sh[f"spd{tag}"]).where(sh[f"spd{tag}"] >= MIN_SPEED_ANGLE)),
        })
    for ax in axs:
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("SADCP speed (m/s)")
        ax.grid(alpha=0.3)
        ax.legend(markerscale=4, fontsize=8)
    axs[0].set_ylabel("|U$_{SADCP}$| − |U$_{geo}$| (m/s)")
    axs[0].set_title("(a) Speed difference vs SADCP speed", loc="left")
    axs[1].set_ylabel("Angle SADCP − geo (°)")
    axs[1].set_ylim(-180, 180)
    axs[1].set_yticks(np.arange(-180, 181, 45))
    axs[1].set_title("(b) Angle difference vs SADCP speed", loc="left")
    fig.tight_layout()
    fig.savefig(DIR_OUT / fname, dpi=200)
    plt.close(fig)
    stats = pd.DataFrame(rows).set_index("SADCP").T
    stats.to_csv(DIR_OUT / "table2_like_stats.csv")
    print(stats.round(3).to_string())
    return stats


# %% Animation — arrows accumulate during a day, cleared at the next day
def make_animation(sh, geo, extent, fname="sadcp_vs_geostrophy_daily.mp4"):
    tag = "100" if ANIM_LEVEL == "layer" else "30"
    lab_ship = (f"SADCP {Z_LAYER[0]:.0f}-{Z_LAYER[1]:.0f} m mean" if tag == "100"
                else f"SADCP {Z_SURF:.0f} m")
    # frames always span the whole cruise, whatever the data coverage
    frames = pd.date_range(pd.Timestamp(CRUISE_START).floor("D"),
                           pd.Timestamp(CRUISE_END), freq=ANIM_STEP)
    print(f"Animation: {len(frames)} frames, {frames[0]:%Y-%m-%d} -> {frames[-1]:%Y-%m-%d %H:%M}")
    a0 = geo.adt.sel(time=frames[0], method="nearest")
    vmin, vmax = np.nanpercentile(geo.adt, [2, 98])

    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.6, 1], hspace=0.35, wspace=0.12)
    ax = fig.add_subplot(gs[:, 0], projection=ccrs.PlateCarree())
    ax_s = fig.add_subplot(gs[0, 1])
    ax_a = fig.add_subplot(gs[1, 1])
    base_map(ax, extent)
    pm = ax.pcolormesh(a0.lon, a0.lat, a0.values, cmap="RdYlBu_r", vmin=vmin, vmax=vmax,
                       shading="auto", transform=ccrs.PlateCarree(), zorder=0)
    plt.colorbar(pm, ax=ax, shrink=0.75, pad=0.02, label="ADT (m)")
    ax.plot(sh.lon, sh.lat, color="0.4", lw=0.5, alpha=0.6, transform=ccrs.PlateCarree(), zorder=2)
    ax.legend(handles=[mlines.Line2D([], [], color="tab:red", lw=2, label=lab_ship),
                       mlines.Line2D([], [], color="tab:blue", lw=2, label="Geostrophic (ADT)"),
                       mlines.Line2D([], [], color="tab:red", lw=1.4, label="AE max-speed contour"),
                       mlines.Line2D([], [], color="tab:blue", lw=1.4, label="CE max-speed contour"),
                       mlines.Line2D([], [], color="k", lw=0.8, ls=":", label="Outermost contour")],
              loc="upper left", fontsize=7, framealpha=0.9).set_zorder(10)

    state = {"day": None, "eddies": [], "cs": None, "dyn": []}
    kw = dict(transform=ccrs.PlateCarree(), scale=ARROW_SCALE, scale_units="inches",
              width=0.003, zorder=6)

    def update(t):
        day = t.floor("D")
        if day != state["day"]:
            state["day"] = day
            a = adt_for_day(geo, day)
            for art in state["eddies"]:
                art.remove()
            state["eddies"] = draw_eddies(ax, day, extent)
            if state["cs"] is not None:
                remove_contour(state["cs"])
                state["cs"] = None
            if a is None:                                   # no ADT map for this day
                pm.set_array(np.full(a0.size, np.nan))
            else:
                pm.set_array(a.values.ravel())
                state["cs"] = ax.contour(a.lon, a.lat, a, levels=np.arange(-2, 3, 0.05),
                                         colors="0.25", linewidths=0.3,
                                         transform=ccrs.PlateCarree(), zorder=1)
        for art in state["dyn"]:
            art.remove()
        state["dyn"] = []

        cur = sh.sel(time=slice(day, t))                      # today only, up to now
        if cur.sizes["time"] > 0:
            state["dyn"] += ax.plot(cur.lon, cur.lat, "k-", lw=1.2,
                                    transform=ccrs.PlateCarree(), zorder=5)
            ok = np.isfinite(cur.lon.values)
            if ok.any():
                i = np.where(ok)[0][-1]
                state["dyn"] += ax.plot(cur.lon.values[i], cur.lat.values[i], marker="*",
                                        ms=13, mfc="yellow", mec="k",
                                        transform=ccrs.PlateCarree(), zorder=8)
            th = cur.resample(time="1h").first()
            X, Y = th.lon.values, th.lat.values
            q1 = ax.quiver(X, Y, th.ug.values, th.vg.values, color="tab:blue", **kw)
            q2 = ax.quiver(X, Y, th[f"u{tag}"].values, th[f"v{tag}"].values,
                           color="tab:red", **kw)
            state["dyn"] += [q1, q2, ax.quiverkey(q2, 0.85, 1.02, 0.5, "0.5 m/s", labelpos="E")]

        # right panels: today's time series
        ax_s.cla()
        ax_a.cla()
        hrs = (pd.to_datetime(cur.time.values) - day) / pd.Timedelta(hours=1)
        ax_s.plot(hrs, cur[f"spd{tag}"], "-", color="tab:red", label=lab_ship)
        ax_s.plot(hrs, cur.spdg, "-", color="tab:blue", label="Geostrophic")
        ax_s.set_xlim(0, 24)
        ax_s.set_ylim(0, max(1.0, float(np.nanmax(sh[f"spd{tag}"])) * 1.05))
        ax_s.set_ylabel("Speed (m/s)")
        ax_s.legend(fontsize=7, loc="upper left")
        ax_s.grid(alpha=0.3)
        da = cur[f"dang{tag}"].where(cur[f"spd{tag}"] >= MIN_SPEED_ANGLE)
        ax_a.plot(hrs, da, "o", ms=3, color="k")
        ax_a.axhline(0, color="k", lw=0.6)
        ax_a.set_xlim(0, 24)
        ax_a.set_ylim(-180, 180)
        ax_a.set_yticks(np.arange(-180, 181, 90))
        ax_a.set_xlabel(f"Hour of {day:%Y-%m-%d} (UTC)")
        ax_a.set_ylabel("Angle SADCP − geo (°)")
        ax_a.grid(alpha=0.3)
        dn = cur[f"dnorm{tag}"]
        ax_s.set_title(f"Today: Δ|U| = {float(dn.mean()):+.2f} ± {float(dn.std()):.2f} m/s"
                       if np.isfinite(dn).any() else "Today: no SADCP yet", fontsize=9)
        ax_a.set_title(f"Δθ = {float(da.mean()):+.0f} ± {float(da.std()):.0f}°"
                       if np.isfinite(da).any() else "", fontsize=9)
        ax.set_title(f"WHIRLS — {t:%Y-%m-%d %H:%M} UTC", loc="left")
        return []

    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000 / ANIM_FPS)
    out = DIR_OUT / fname

    def progress(i, n):
        if i % 12 == 0 or i == n - 1:
            print(f"  frame {i + 1}/{n}  ({frames[i]:%Y-%m-%d %H:%M})", flush=True)

    if animation.writers.is_available("ffmpeg"):
        anim.save(out, writer=animation.FFMpegWriter(fps=ANIM_FPS, bitrate=3000), dpi=120,
                  progress_callback=progress)
    else:
        out = out.with_suffix(".gif")
        anim.save(out, writer=animation.PillowWriter(fps=ANIM_FPS), dpi=90,
                  progress_callback=progress)
    plt.close(fig)
    print(f"  saved {out}")
    return out


# %% Run figures
if __name__ == "__main__" or "get_ipython" in globals():
    prof = profile_stats(sh)
    if dp is not None:
        prof = pd.concat([prof.assign(instrument="OS150"),
                          profile_stats(dp).assign(instrument="OS38")])
    prof.to_csv(DIR_OUT / "stats_vs_depth.csv", index=False)

    plot_fig7(sh, dp, geo, CRUISE_START, CRUISE_END, "fig7_like_whole_cruise.png",
              title="Whole cruise", extent=EXTENT)
    for name, (a, b) in SECTIONS.items():
        plot_fig7(sh, dp, geo, a, b, f"fig7_like_{name}.png", title=name)
    plot_fig8_and_stats(sh)
    make_animation(sh, geo, EXTENT)