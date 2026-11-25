import xarray as xr
import numpy as np
import netCDF4 as nc
import xarray as xr
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import numpy as np
import pandas as pd
import gsw
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import cmocean.cm as cmo
import matplotlib.dates as mdates
import matplotlib.lines as mlines
from matplotlib.colors import TwoSlopeNorm
from pyproj import Geod
import warnings
import cmocean
from matplotlib.path import Path as MplPath
import re
 
###### CONSTANTS ######
SADCP_AVG = "10min"          # SADCP time averaging (~3 km at 10 kn)
USE_DETIDED = True
Z_SPLIT = 150.0              # m: OS150 above, OS38 below
ZMAX_SADCP = 1000.0
Z_FILL_GAP = 40.0            # m: fill vertical gaps up to this size (seam between OS150/OS38)
DS_KM = 3.0                  # bin size along the section line (panel c)
HEADING_SMOOTH_KM = 2.0      # track length used for the local heading (panel b)
CENTER_SMOOTH_DAYS = 3       # running mean on the daily TOEddies centers
R_KM = 6371.0
G = 9.81
RHO0 = 1025.0
OMEGA = 7.2921e-5
geod = Geod(ellps="WGS84")
 
ADT_MIN, ADT_MAX = -0.3, 1.2
adt_norm = TwoSlopeNorm(vmin=ADT_MIN, vcenter=0, vmax=ADT_MAX)
levels = np.arange(ADT_MIN, ADT_MAX + 0.001, 0.02)
adt_ticks = [-0.3, 0, 0.3, 0.6, 0.9, 1.2]
NEN_COLOR = "tab:green"
D_DEPTHS = (50, 100, 200, 400, 700)   # m: depths drawn on panel (d) (nearest SADCP level)
C_LINE_STEP = 0.1                      # m/s: interval of the black contour lines on panel (c)
N_FILL = 50    
NEN_SEARCH_KM = 80.0         # candidate centers: square of +- this around the TOEddies center
NEN_GRID = 81                # candidate grid points per side (81 -> 2 km spacing)
NEN_MIN_FRAC = 0.5           # skip depth levels with fewer valid fixes than this fraction
NEN_LAYER = (30.0, 200.0)    # m: layer for the single robust center (map + print)
CENTER_BY_DEPTH = True       # True: green line = center at every depth; False: layer center only
SECTION_X_FROM = "start" # m: layer for the single robust center (map + print)
CRUISE_START = "2026-05-14"
CRUISE_END   = "2026-07-27"
 
###### PATHS ######
BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
DIR_SADCP = BASE / "data/SADCP"
SADCP_FILES = {"OS150": DIR_SADCP / "WHIRLS_NB150_TOT_osite.nc",   # full cruise (README)
               "OS38": DIR_SADCP / "WHIRLS_BB38_TOT_osite.nc"}
DIR_LADCP = BASE / "data/LADCP/traitement_avec_SADCP"
DIR_SSH = BASE / "data/SSH"
DIR_EDDIES = {"TOEDDIES_better_ADT": BASE / "data/TOEDDIES/TOEDDIES_better_ADT",
              "TOEDDIES": BASE / "data/TOEDDIES/TOEDDIES"}
FNAME_BATHY = BASE / "data/bathymetry/etopo2.nc"
DIR_TRACKS = BASE / "figures/eddy_tracks"
DIR_CODE = BASE / "analysis"   
 
DIR_DRIFTERS =  BASE / "data/DRIFTERS/drifters"
DIR_GEOJSON = BASE / "data/DRIFTERS/drifters_tracks.geojson"
 
DIR_MVP = BASE / "data/MVP_Data_processed"
 
CMAPS = dict(CT=cmo.thermal, SA=cmo.haline, sig=cmo.dense, pv=cmo.balance)
COLS = dict(p="PRES (dbar)", t="TEMP (°C)", sp="SAL (psu)", c="COND (mS/cm)",
            time="TIME", lat="LAT (°N)", lon="LON (°E)")
# Bathymetry
lon_min, lon_max, lat_min, lat_max = 8, 25, -41, -34
extent = [lon_min, lon_max, lat_min, lat_max]
bathymetry_data = xr.open_dataset(FNAME_BATHY)
lon_subset = bathymetry_data["lon"].sel(lon=slice(lon_min, lon_max))
lat_subset = bathymetry_data["lat"].sel(lat=slice(lat_min, lat_max))
depth_subset = bathymetry_data["topo"].sel(lon=slice(lon_min, lon_max),
                                           lat=slice(lat_min, lat_max))
 
# File lists
daily_ssh_files = sorted(os.listdir(DIR_SSH))
eddy_files_better = sorted(os.listdir(DIR_EDDIES['TOEDDIES_better_ADT']))
eddy_files_std = sorted(os.listdir(DIR_EDDIES['TOEDDIES']))
 
###### SADCP FUNCTIONS ######
def to_xy(lon, lat, lon0, lat0):
    """Local tangent plane, km."""
    x = np.deg2rad(np.asarray(lon, float) - lon0) * R_KM * np.cos(np.deg2rad(lat0))
    y = np.deg2rad(np.asarray(lat, float) - lat0) * R_KM
    return x, y
 
 
def to_lonlat(x, y, lon0, lat0):
    lon = lon0 + np.rad2deg(np.asarray(x, float) / (R_KM * np.cos(np.deg2rad(lat0))))
    lat = lat0 + np.rad2deg(np.asarray(y, float) / R_KM)
    return lon, lat
 
def to_sn(lo, la, lon0, lat0, origin, es, en):
        px, py = to_xy(lo, la, lon0, lat0)
        return ((px - origin[0]) * es[0] + (py - origin[1]) * es[1],
                (px - origin[0]) * en[0] + (py - origin[1]) * en[1])
 
def _find(ds, names, required=True):
    for n in names:
        if n in ds.variables:
            return n
    low = {k.lower(): k for k in ds.variables}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    if required:
        raise KeyError(f"None of {names} in {list(ds.variables)}")
    return None
 
def cascade_add_time(sadcp):
    """Turn JULD (days since REFERENCE_DATE_TIME) into a datetime 'time' coordinate."""
    juld = sadcp["JULD"]
    tdim = juld.dims[0]
    if np.issubdtype(juld.dtype, np.datetime64):
        t = juld.values
    else:
        ref = sadcp["REFERENCE_DATE_TIME"].values
        ref = ref.tobytes().decode() if ref.dtype.kind in "SU" or (ref.dtype == object and ref.ndim) else str(ref)
        ref = "".join(ch for ch in ref if ch.isdigit()) or "19500101000000"
        ref = pd.to_datetime(ref[:14], format="%Y%m%d%H%M%S")
        t = ref + pd.to_timedelta(juld.values, unit="D")
    sadcp = sadcp.assign_coords(time=(tdim, pd.DatetimeIndex(t).values))
    sadcp = sadcp.swap_dims({tdim: "time"})
 
    # attach bin depths as a 1-D 'depth' coordinate on the depth-bin dimension
    zdim = [d for d in sadcp["UVEL_ADCP"].dims if d != "time"][0]
    dep = sadcp["DEPH"]
    if "time" in dep.dims:                       # time x depth -> one vector
        dep = dep.median("time", skipna=True)
    sadcp = sadcp.assign_coords(depth=(zdim, np.abs(dep.values)))
    return sadcp
 
 
def sadcp_for_transect(data, sadcp, pad="5min"):
    """SADCP subset covering the MVP transect's time window (+/- pad)."""
    t = pd.to_datetime(data["TIME"])
    pad = pd.Timedelta(pad)
    return sadcp.sel(time=slice(t.min() - pad, t.max() + pad))
 
 
def sadcp_on_mvp_times(data, sadcp_tr, max_gap="10min",
                       keep=("UVEL_ADCP", "VVEL_ADCP", "UVEL_ADCP_CORTIDE",
                             "VVEL_ADCP_CORTIDE", "LATITUDE", "LONGITUDE", "PGOOD_ADCP")):
    """Interpolate selected SADCP variables onto each MVP sample time (keeps 'depth' coord).
    MVP samples more than `max_gap` from an ADCP ensemble are set to NaN."""
    keep = [v for v in keep
            if v in sadcp_tr and "time" in sadcp_tr[v].dims
            and np.issubdtype(sadcp_tr[v].dtype, np.number)]
    ds = sadcp_tr[keep]                          # 'depth' coordinate comes along
 
    t_mvp = pd.to_datetime(data["TIME"]).values
    out = ds.interp(time=t_mvp)
 
    t_adcp = ds.time.values
    idx = np.clip(np.searchsorted(t_adcp, t_mvp), 1, len(t_adcp) - 1)
    gap = np.minimum(np.abs(t_mvp - t_adcp[idx - 1]), np.abs(t_adcp[idx] - t_mvp))
    ok = xr.DataArray(gap <= pd.Timedelta(max_gap), dims="time", coords={"time": t_mvp})
    return out.where(ok)
 
 
def load_sadcp_raw(path):
    """Full-resolution SADCP as Dataset(time, z) with u, v, lat, lon (QC flag == 1)."""
    ds = xr.open_dataset(path)
    tn = _find(ds, ["JULD_ADCP", "TIME", "time", "JULD"])
    zn = _find(ds, ["DEPH", "DEPTH", "depth"])
    un = _find(ds, ["UVEL_ADCP_CORTIDE"], required=False) if USE_DETIDED else None
    vn = _find(ds, ["VVEL_ADCP_CORTIDE"], required=False) if USE_DETIDED else None
    if un is None or vn is None:
        un, vn = _find(ds, ["UVEL_ADCP"]), _find(ds, ["VVEL_ADCP"])
    fn = _find(ds, ["CAS_CURRENT_FLAG"], required=False)
    print(f"  {Path(path).name}: using {un}, {vn}")
    tdim = ds[tn].dims[0]
    U = ds[un].squeeze(drop=True)
    zdim = [d for d in U.dims if d != tdim][0]
    u = U.transpose(tdim, zdim).values.astype(float)
    v = ds[vn].squeeze(drop=True).transpose(tdim, zdim).values.astype(float)
    traw = ds[tn].values
    time = (pd.to_datetime(traw) if np.issubdtype(traw.dtype, np.datetime64)
            else pd.Timestamp("1950-01-01") + pd.to_timedelta(traw.astype(float), unit="D"))
    if fn is not None:
        fl = ds[fn].squeeze(drop=True)
        fl = fl.transpose(tdim, zdim).values if zdim in fl.dims else fl.values[:, None]
        good = (fl == 1) if np.issubdtype(fl.dtype, np.number) else \
            np.char.strip(np.char.decode(fl.astype("S"), "ascii", "ignore")) == "1"
        good = np.broadcast_to(good, u.shape)
        u, v = np.where(good, u, np.nan), np.where(good, v, np.nan)
    u[np.abs(u) > 5] = np.nan
    v[np.abs(v) > 5] = np.nan
    zda = ds[zn].squeeze(drop=True)
    z = np.abs(np.nanmedian(zda.transpose(tdim, zdim).values, axis=0) if zda.ndim == 2
               else zda.values).astype(float)
    keepz = np.isfinite(z)
    lat = ds[_find(ds, ["LATITUDE", "lat"])].squeeze(drop=True).values.astype(float)
    lon = ds[_find(ds, ["LONGITUDE", "lon"])].squeeze(drop=True).values.astype(float)
    lat = np.where(np.abs(lat) <= 90, lat, np.nan)
    lon = np.where(np.abs(lon) <= 360, lon, np.nan)
    lon = np.where(lon > 180, lon - 360, lon)
    out = xr.Dataset({"u": (("time", "z"), u[:, keepz]), "v": (("time", "z"), v[:, keepz]),
                      "lat": ("time", lat), "lon": ("time", lon)},
                     coords={"time": np.asarray(time, "datetime64[ns]"), "z": z[keepz]})
    out = out.isel(time=np.asarray(pd.notna(time))).sortby("time").sortby("z")
    out = out.isel(time=~pd.Index(out.time.values).duplicated())
    ds.close()
    return out
 
def sadcp_window(inst, sadcp=None, t0=None, t1=None, files=None):
    """OS150 above Z_SPLIT + OS38 below, averaged over SADCP_AVG.
    t0, t1 : optional time limits; None = from the start / to the end of the record.
    files  : dict {"OS150": path, "OS38": path}; default SADCP_FILES.
    Returns a dict of numpy arrays: time, lon, lat, z, u (nt, nz), v (nt, nz)."""
    files = SADCP_FILES if files is None else files
    parts = []
    if sadcp is None:
        for inst, f in files.items():
            if not Path(f).exists():
                print(f"  ! missing {f}")
                continue
            d = load_sadcp_raw(f).sel(time=slice(t0, t1)).resample(time=SADCP_AVG).mean()
            zz = d.z.values
            keep = zz <= Z_SPLIT if inst == "OS150" else (zz > Z_SPLIT) & (zz <= ZMAX_SADCP)
            parts.append(d.isel(z=np.where(keep)[0]))
        if not parts:
            raise FileNotFoundError("no SADCP file found")
    else:
        print("Calculating Mean")
        sadcp = sadcp.resample(time=SADCP_AVG).mean() 
        zz = sadcp.z.values
        keep = zz <= Z_SPLIT if inst == "OS150" else (zz > Z_SPLIT) & (zz <= ZMAX_SADCP)
        parts = [sadcp.isel(z=np.where(keep)[0])]
    out = parts[0]
    for p in parts[1:]:
        # outer join: keep times where either instrument has data
        out, p = xr.align(out, p, join="outer", exclude=["z"])
        out = xr.Dataset({"u": xr.concat([out.u, p.u], "z"), "v": xr.concat([out.v, p.v], "z"),
                          "lat": out.lat.fillna(p.lat), "lon": out.lon.fillna(p.lon)})
    out = out.sortby("z")
    ok = np.isfinite(out.lat.values) & np.isfinite(out.lon.values) & np.isfinite(out.u.values).any(1)
    out = out.isel(time=ok)
    return dict(time=out.time.values, lon=out.lon.values, lat=out.lat.values, z=out.z.values,
                u=out.u.values, v=out.v.values)
 
 
 
def track_heading(lon, lat, smooth_km=HEADING_SMOOTH_KM):
    """Ship course (deg clockwise from north) over +- smooth_km/2 of track; NaN on station."""
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    x = along_track_distance(lon, lat) / 1000.0
    ok = np.isfinite(x) & np.isfinite(lon) & np.isfinite(lat)
    heading = np.full(lon.shape, np.nan)
    if ok.sum() < 2:
        return heading
    xs, first = np.unique(x[ok], return_index=True)
    lons, lats = lon[ok][first], lat[ok][first]
    if xs.size < 2:
        return heading
    half = smooth_km / 2.0
    xa = np.clip(x[ok] - half, xs[0], xs[-1])
    xb = np.clip(x[ok] + half, xs[0], xs[-1])
    az, _, d = geod.inv(np.interp(xa, xs, lons), np.interp(xa, xs, lats),
                        np.interp(xb, xs, lons), np.interp(xb, xs, lats))
    heading[ok] = np.where(d > 250.0 * smooth_km, az, np.nan) % 360.0
    return heading
 
def fit_line(x, y):
    """Principal axis of the points: origin, along unit vector (direction of travel),
    normal unit vector (LEFT of travel)."""
    ok = np.isfinite(x) & np.isfinite(y)
    xm, ym = x[ok].mean(), y[ok].mean()
    w, vec = np.linalg.eigh(np.cov(np.vstack([x[ok] - xm, y[ok] - ym])))
    es = vec[:, np.argmax(w)]
    xi, yi = x[ok][[0, -1]], y[ok][[0, -1]]
    if (xi[1] - xi[0]) * es[0] + (yi[1] - yi[0]) * es[1] < 0:
        es = -es
    return np.array([xm, ym]), es, np.array([-es[1], es[0]])
 
def along_track_distance(lon, lat):
    """Cumulative geodesic distance along the track [m]."""
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    ok = np.isfinite(lon) & np.isfinite(lat)
    dist = np.full(lon.shape, np.nan)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return dist
    _, _, seg = geod.inv(lon[idx[:-1]], lat[idx[:-1]], lon[idx[1:]], lat[idx[1:]])
    dist[idx] = np.concatenate([[0.0], np.cumsum(seg)])
    return dist
 
def orthogonal_velocity(lon, lat, u, v):
    """Velocity perpendicular to the ship's local course, + to the LEFT of travel."""
    az = np.deg2rad(track_heading(lon, lat))
    nx, ny = -np.cos(az)[:, None], np.sin(az)[:, None]
    return u * nx + v * ny
 
def bin_columns(coord, vals, edges, nmin=1):
    """Mean of vals (n, nz) in bins of coord. Returns centers, (nbins, nz) mean, counts."""
    idx = np.digitize(coord, edges) - 1
    out = np.full((len(edges) - 1, vals.shape[1]), np.nan)
    cnt = np.zeros(len(edges) - 1, int)
    for b in range(len(edges) - 1):
        m = idx == b
        cnt[b] = m.sum()
        if cnt[b] >= nmin:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out[b] = np.nanmean(vals[m], axis=0)
    return 0.5 * (edges[:-1] + edges[1:]), out, cnt
 
 
def fill_z_gaps(V, z, max_gap=Z_FILL_GAP):
    """Linear interpolation in depth across interior gaps no thicker than max_gap."""
    V = V.copy()
    for i in range(V.shape[0]):
        ok = np.isfinite(V[i])
        if ok.sum() < 2:
            continue
        zk = z[ok]
        bad = np.where(~ok & (z > zk.min()) & (z < zk.max()))[0]
        if not len(bad):
            continue
        j = np.searchsorted(zk, z[bad])
        small = (zk[j] - zk[j - 1]) <= max_gap
        V[i, bad[small]] = np.interp(z[bad[small]], zk, V[i, ok])
    return V
 
def cascade_to_window_format(ds, detided=True, pgood_min=None):
    """
    Convert a CASCADE-layout SADCP dataset (e.g. sadcp_tr or sadcp_mvp, after
    cascade_add_time) to the layout sadcp_window expects:
        Dataset(time, z) with variables u, v (time, z) and lat, lon (time),
        z = positive depth in m, sorted increasing.
    The input dataset is not modified, so the rest of the notebook keeps working.
 
    detided   : use UVEL/VVEL_ADCP_CORTIDE if present, else UVEL/VVEL_ADCP.
    pgood_min : if set (e.g. 50), mask u, v where PGOOD_ADCP < pgood_min.
    """
    if detided and "UVEL_ADCP_CORTIDE" in ds and "VVEL_ADCP_CORTIDE" in ds:
        un, vn = "UVEL_ADCP_CORTIDE", "VVEL_ADCP_CORTIDE"
    else:
        un, vn = "UVEL_ADCP", "VVEL_ADCP"
 
    zdim = [d for d in ds[un].dims if d != "time"][0]          # e.g. "N_LEVEL"
 
    # depth of each bin: 'depth' coord from cascade_add_time, else DEPH
    if "depth" in ds.coords:
        z = ds["depth"].values
    else:
        dep = ds["DEPH"]
        z = dep.median("time", skipna=True).values if "time" in dep.dims else dep.values
    z = np.abs(np.asarray(z, dtype=float))
 
    u = ds[un].transpose("time", zdim).values.astype(float)
    v = ds[vn].transpose("time", zdim).values.astype(float)
    if pgood_min is not None and "PGOOD_ADCP" in ds:
        pg = ds["PGOOD_ADCP"].transpose("time", zdim).values
        if np.nanmax(pg) > 0:                      # OS38 file: PGOOD not filled (all 0) -> skip
            u[pg < pgood_min] = np.nan
            v[pg < pgood_min] = np.nan
        else:
            print("  ! PGOOD_ADCP is all 0 in this file -> PGOOD mask skipped")
    # same sanity limit as load_sadcp_raw: drop fill values / spikes
    u[np.abs(u) > 5] = np.nan
    v[np.abs(v) > 5] = np.nan
 
    lon = ds["LONGITUDE"].values.astype(float)
    lon = np.where(lon > 180, lon - 360, lon)                  # same convention as load_sadcp_raw
 
    out = xr.Dataset(
        {"u": (("time", "z"), u), "v": (("time", "z"), v),
         "lat": ("time", ds["LATITUDE"].values.astype(float)), "lon": ("time", lon)},
        coords={"time": ds["time"].values, "z": z})
    out = out.isel(z=np.where(np.isfinite(z))[0]).sortby("z")
    return out
 
# =============================================================================
# Nencioli et al. (2008): your algorithm, corrected and vectorised
# =============================================================================
def nencioli_center(lon, lat, u, v, lon_c0, lat_c0, search_km=NEN_SEARCH_KM, n_grid=NEN_GRID,
                    min_frac=NEN_MIN_FRAC):
    """
    Eddy center from the velocities measured along a ship track (Nencioli et al. 2008).
 
    For every candidate center (lon2, lat2) of a square grid (your WL..EL, SL..NL box,
    here +- search_km around lon_c0, lat_c0), every measured velocity is split into
    radial and tangential components about that candidate. The center is the candidate
    where the RMS of the radial velocity is smallest. Done for each depth at once.
 
    lon, lat : (nfix,)          positions of the measurements
    u, v     : (nfix, nz)       velocities (east, north)
    Returns a dict of (nz,) arrays: lon, lat (center), rms (min radial RMS), mtg (mean
    tangential velocity there, + anticlockwise), edge (min on the grid edge), ok (enough
    data), plus the grid and the full RMS field rms_grid (n_grid, n_grid, nz).
 
    Fixes compared with the original loop:
      * angles in radians (cos/sin were given degrees);
      * dx, dy are signed east/north distances from the candidate to the measurement,
        so arctan2 gives the right angle in every quadrant (no quadrant corrections);
      * radial, tangent, Elliptic_RMS, Elliptic_Mtg are separate arrays (they were the
        same array: `radial = A; tangent = radial` makes them aliases);
      * the minimum is located with unravel_index on the 2-D grid (np.where(...)[0] only
        gave the row);
      * the RMS is divided by the number of valid values, not the full length;
      * the SADCP arrays are not modified in place.
    """
    R, f = R_KM, np.pi / 180
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    u, v = np.atleast_2d(np.asarray(u, float).T).T, np.atleast_2d(np.asarray(v, float).T).T
 
    # candidate grid
    dlat = search_km / (R * f)
    dlon = dlat / np.cos(lat_c0 * f)
    lon2 = np.linspace(lon_c0 - dlon, lon_c0 + dlon, n_grid)
    lat2 = np.linspace(lat_c0 - dlat, lat_c0 + dlat, n_grid)
    LON, LAT = np.meshgrid(lon2, lat2)                       # (ny, nx)
    lc, la = LON.ravel()[:, None], LAT.ravel()[:, None]      # (ncand, 1)
 
    # signed distances candidate -> measurement [km] and angle from east, anticlockwise
    dx = (lon[None, :] - lc) * f * R * np.cos(la * f)
    dy = (lat[None, :] - la) * f * R
    angle = np.arctan2(dy, dx)                               # (ncand, nfix)
    C, S = np.cos(angle), np.sin(angle)
 
    # radial = U cos + V sin, tangential = -U sin + V cos, for every candidate and depth
    valid = np.isfinite(u) & np.isfinite(v) & (np.isfinite(lon) & np.isfinite(lat))[:, None]
    U, V = np.where(valid, u, 0.0), np.where(valid, v, 0.0)
    C, S = np.nan_to_num(C), np.nan_to_num(S)
    cnt = valid.sum(axis=0).astype(float)                    # valid values per depth
    with np.errstate(invalid="ignore", divide="ignore"):
        rad2 = (C ** 2) @ (U ** 2) + 2 * (C * S) @ (U * V) + (S ** 2) @ (V ** 2)
        rms = np.sqrt(rad2 / cnt)                            # (ncand, nz)
        mtg = ((-S) @ U + C @ V) / cnt
 
    nz = u.shape[1]
    ok = cnt >= min_frac * np.isfinite(lon).sum()
    k = np.zeros(nz, int)
    k[ok] = np.nanargmin(rms[:, ok], axis=0)
    iy, ix = np.unravel_index(k, LON.shape)
    edge = (iy == 0) | (iy == n_grid - 1) | (ix == 0) | (ix == n_grid - 1)
    nan = lambda a: np.where(ok, a, np.nan)
    return dict(lon=nan(lon2[ix]), lat=nan(lat2[iy]), rms=nan(rms[k, np.arange(nz)]),
                mtg=nan(mtg[k, np.arange(nz)]), edge=edge & ok, ok=ok,
                grid_lon=lon2, grid_lat=lat2, rms_grid=rms.reshape(n_grid, n_grid, nz))
 
###### TOEDDIES FUNCTIONS ######
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
 
def get_eddy_files(date_compact, date):
    """Use the better-ADT TOEddies when available for that day, else the standard ones."""
    files = [DIR_EDDIES['TOEDDIES_better_ADT'] / f for f in eddy_files_better if date_compact in f]
    if files:
        return files, "TOEddies (better ADT)"
    files = [DIR_EDDIES['TOEDDIES'] / f for f in eddy_files_std if date in f]
    if files:
        return files, "TOEddies (standard)"
    return [], None
 
 
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
from shapely.geometry import LineString
 
def eddy_on_transect(eddy_files, lon_tr, lat_tr):
    """Pick the TOEddies eddy whose OUTER contour the transect crosses
    (nearest centre if several). Returns dict(center, maxspeed, outer, cyc) or None."""
    if not eddy_files:
        return None
    anti, cyc = load_eddies(eddy_files)
    track = LineString(np.column_stack([lon_tr, lat_tr]))
    lat0 = np.nanmean(lat_tr)
    best = None
    for eddies, is_cyc in [(anti, False), (cyc, True)]:
        ms, oc, cent = eddies["maxspeed"], eddies["outcontour"], eddies["centers"]
        for i in range(len(cent)):
            if ms.geometry[i] is None or oc.geometry[i] is None:
                continue
            if not track.intersects(oc.geometry[i]):          # transect must cross the eddy
                continue
            cx, cy = cent.geometry[i].x, cent.geometry[i].y
            # distance (km) from the centre to the nearest point of the transect
            d = np.nanmin(np.hypot((lon_tr - cx) * 111.32 * np.cos(np.deg2rad(lat0)),
                                   (lat_tr - cy) * 111.32))
            if best is None or d < best[0]:
                best = (d, dict(center=(cx, cy),
                                maxspeed=np.array(ms.geometry[i].exterior.xy).T,  # (N, 2) lon, lat
                                outer=np.array(oc.geometry[i].exterior.xy).T,
                                cyc=is_cyc))
    return None if best is None else best[1]
 
 
def eddies_on_transect(eddy_files, lon_tr, lat_tr):
    """ALL TOEddies eddies whose OUTER contour the transect crosses (cyclones and
    anticyclones), ordered along the track (first crossed first).
    Returns a list of dict(center, maxspeed, outer, cyc, dist_km); [] if none."""
    if not eddy_files:
        return []
    lon_tr, lat_tr = np.asarray(lon_tr, float), np.asarray(lat_tr, float)
    ok = np.isfinite(lon_tr) & np.isfinite(lat_tr)
    lon_tr, lat_tr = lon_tr[ok], lat_tr[ok]
    anti, cyc = load_eddies(eddy_files)
    track = LineString(np.column_stack([lon_tr, lat_tr]))
    lat0 = np.nanmean(lat_tr)
    found = []
    for eddies, is_cyc in [(anti, False), (cyc, True)]:
        ms, oc, cent = eddies["maxspeed"], eddies["outcontour"], eddies["centers"]
        for i in range(len(cent)):
            if ms.geometry[i] is None or oc.geometry[i] is None:
                continue
            if not track.intersects(oc.geometry[i]):          # transect must cross the eddy
                continue
            cx, cy = cent.geometry[i].x, cent.geometry[i].y
            d = np.hypot((lon_tr - cx) * 111.32 * np.cos(np.deg2rad(lat0)), (lat_tr - cy) * 111.32)
            found.append(dict(center=(cx, cy),
                              maxspeed=np.array(ms.geometry[i].exterior.xy).T,  # (N, 2) lon, lat
                              outer=np.array(oc.geometry[i].exterior.xy).T,
                              cyc=is_cyc,
                              dist_km=float(np.nanmin(d)),       # centre -> closest track point
                              _order=int(np.nanargmin(d))))      # where along the track
    found.sort(key=lambda e: e.pop("_order"))
    return found
 
 
def crossings_along_track(lon, lat, x, xy):
    """x where the ship track enters / leaves contour xy (midpoint of the two fixes)."""
    ok = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(x)
    if xy is None or ok.sum() < 2:
        return np.array([])
    idx = np.flatnonzero(ok)
    io = MplPath(xy).contains_points(np.c_[lon[idx], lat[idx]])
    k = np.flatnonzero(io[1:] != io[:-1])
    return 0.5 * (x[idx[k]] + x[idx[k + 1]])
 
 
def crossings_on_line(xy, to_sn):
    """s where contour xy crosses the straight section line (n = 0)."""
    if xy is None:
        return np.array([])
    s, n = to_sn(xy[:, 0], xy[:, 1])
    k = np.flatnonzero(np.sign(n[:-1]) * np.sign(n[1:]) < 0)
    return s[k] - n[k] * (s[k + 1] - s[k]) / (n[k + 1] - n[k])
 
 
 
###### MVP FUNCTIONS ######
 
def make_to_sn(sec):
    """(lon, lat) -> (s, n) in km on the MVP section axis:
    s = along-track from the first cast (same as sec["x"]/1e3), n = to the left."""
    lon0, lat0, ex = sec["lon"][0], sec["lat"][0], sec["xdir"]
    R = 6371e3
    def to_sn(lo, la):
        lo, la = np.asarray(lo, float), np.asarray(la, float)
        e = np.deg2rad(lo - lon0) * R * np.cos(np.deg2rad(lat0))
        nn = np.deg2rad(la - lat0) * R
        return (e * ex[0] + nn * ex[1]) / 1e3, (-e * ex[1] + nn * ex[0]) / 1e3
    return to_sn
 
# ---------------------------------------------------------------------------
# One call per transect
# ---------------------------------------------------------------------------
def subset_section(sec, keep):
    """Keep only the casts where `keep` is True (every array whose last axis is the casts)."""
    keep = np.asarray(keep, bool)
    m = keep.size
    out = {}
    for k, val in sec.items():
        if isinstance(val, pd.DatetimeIndex):
            out[k] = val[keep]
        elif isinstance(val, np.ndarray) and val.ndim >= 1 and val.shape[-1] == m:
            out[k] = val[..., keep]
        else:
            out[k] = val
    return out


def auto_flip(xdir):
    """True if the x axis must be inverted so the section reads like the map:
    zonal transect -> West on the left, East on the right;
    meridional transect -> South on the left, North on the right.
    xdir = (east, north) unit vector of +x (= direction the ship moves)."""
    east, north = xdir
    return bool(east < 0) if abs(east) >= abs(north) else bool(north < 0)


def process_transect(df, nav=None, adcp=None, dp=1.0, smooth_z=5, smooth_x=3,
                     min_cast_depth=20.0, cast_col=None, title="", savepath=None,
                     map_draw=None, map_args=None, eddy=None, x_max_km=None, x_min_km=None,
                     **plot_kw):
    """
    df   : raw MVP DataFrame for one transect
    nav  : (time, lon, lat) of ship GPS, used if MVP LAT/LON are NaN
    cast_col : column identifying casts (e.g. "CAST" = source file). If None,
           casts are detected from the pressure record (needs one continuous time series).
    map_draw : function(ax, crs) drawing the map background. If given, the figure is
           map + Theta/S_A/sigma0; if None, it is Theta/S_A/sigma0/PV (PV is
           computed either way: sec["q"]; plot it with plot_section(sec)).
    adcp : (lon, lat, depth, u, v) for zeta, e.g.
           (adcp_cut.lon_adcp.values, adcp_cut.lat_adcp.values, adcp_cut.depth.values,
            adcp_cut.uvel_adcp.values, adcp_cut.vvel_adcp.values)
    x_min_km, x_max_km : keep only the casts between these along-track distances (km from
           the first cast, in the direction of motion), e.g. x_max_km=105 to drop
           everything after a gap. x is then recomputed on the kept casts.
    flip_x (in plot_kw): "auto" (default with a map) -> the section reads like the map
           (W left / E right, or S left / N right); True / False to force it.
    """
    if nav is not None:
        df = fill_positions(df, *nav)
    if df[COLS["lat"]].isna().all():
        raise ValueError("No positions: pass nav=(time, lon, lat) from the ship GPS / ADCP.")
    if cast_col is not None:
        casts = split_by_cast_id(df, cast_col)
    else:
        casts = split_downcasts(df, min_cast_depth=min_cast_depth)
    print(f"{len(casts)} casts used")
    casts = fill_cast_positions(casts)
    sec = build_section(casts, dp=dp)
    sec["x"], sec["xdir"] = along_track_x(sec["lon"], sec["lat"])

    if x_min_km is not None or x_max_km is not None:
        xkm = sec["x"] / 1e3
        keep = np.ones(xkm.size, bool)
        if x_min_km is not None:
            keep &= xkm >= x_min_km
        if x_max_km is not None:
            keep &= xkm <= x_max_km
        print(f"cut {x_min_km}-{x_max_km} km: {keep.sum()} of {keep.size} casts kept")
        sec = subset_section(sec, keep)
        # origin = first kept cast, direction = first -> last kept cast
        sec["x"], sec["xdir"] = along_track_x(sec["lon"], sec["lat"])
 
    if map_draw is not None:
        fig, axs = plot_section_with_map(sec, adcp, map_draw=map_draw, map_args=map_args, title=title,
                                         savepath=savepath, edd=eddy, **plot_kw)
    else:
        fig, axs = plot_section(sec, title=title, savepath=savepath, **plot_kw)
    return sec, fig, axs
 
def fill_positions(df, nav_time, nav_lon, nav_lat):
    """Interpolate ship lon/lat (e.g. new_adcp.time, lon_adcp, lat_adcp) onto MVP times."""
    df = df.copy()
    t_mvp = pd.to_datetime(df[COLS["time"]]).values.astype("datetime64[ns]").astype("int64")
    t_nav = np.asarray(nav_time, dtype="datetime64[ns]").astype("int64")
    o = np.argsort(t_nav)
    for col, val in [(COLS["lon"], nav_lon), (COLS["lat"], nav_lat)]:
        filled = np.interp(t_mvp, t_nav[o], np.asarray(val, float)[o], left=np.nan, right=np.nan)
        df[col] = df[col].fillna(pd.Series(filled, index=df.index))
    return df
 
def check_units(df):
    """SAL should match Practical Salinity recomputed from conductivity."""
    sp_from_c = gsw.SP_from_C(df[COLS["c"]], df[COLS["t"]], df[COLS["p"]])
    diff = (df[COLS["sp"]] - sp_from_c).abs()
    print(f"|SAL - SP_from_C| median = {diff.median():.4f}, 95th pct = {diff.quantile(.95):.4f}")
    print("-> < ~0.01 means SAL is Practical Salinity and COND is in mS/cm.")
 
def split_downcasts(df, min_cast_depth=20.0, smooth_n=11):
    """Return a list of DataFrames, one per downcast (surface -> bottom of fish drop)."""
    df = df.sort_values(COLS["time"]).reset_index(drop=True)
    p = df[COLS["p"]].rolling(smooth_n, center=True, min_periods=1).median().values
    bottoms, _ = find_peaks(p, prominence=min_cast_depth)
    tops, _ = find_peaks(-p, prominence=min_cast_depth)
    casts = []
    for b in bottoms:
        prev_tops = tops[tops < b]
        start = prev_tops[-1] if len(prev_tops) else 0
        cast = df.iloc[start:b + 1]
        # keep only samples where the fish is actually descending
        dp = np.gradient(p[start:b + 1])
        cast = cast[dp > 0]
        if len(cast) > 10:
            casts.append(cast)
    return casts
 
def split_by_cast_id(df, cast_col="CAST", downcast_only=True):
    """
    One cast per value of `cast_col` (e.g. one processed MVP CSV file = one cast).
    downcast_only: within each cast keep samples from the start down to max pressure.
    """
    casts = []
    for _, c in df.groupby(cast_col, sort=False):
        c = c.sort_values(COLS["time"])
        if downcast_only:
            imax = np.argmax(c[COLS["p"]].values)
            if imax > 10:
                c = c.iloc[:imax + 1]
        if len(c) > 10:
            casts.append(c)
    # order casts in time so x increases along the transect
    casts.sort(key=lambda c: pd.to_datetime(c[COLS["time"]]).iloc[0])
    return casts
 
def fill_cast_positions(casts):
    """
    Casts whose LAT/LON are all NaN (T/S fine, GPS missing) get a position
    interpolated in time from the neighbouring casts (linear extrapolation at the
    ends of the transect). Rows with NaN LAT/LON inside a cast need nothing:
    the cast position is the mean of the valid rows.
    """
    t = np.array([pd.to_datetime(c[COLS["time"]]).iloc[0].value for c in casts], float)
    lat = np.array([c[COLS["lat"]].mean() for c in casts])
    lon = np.array([c[COLS["lon"]].mean() for c in casts])
    ok = np.isfinite(lat) & np.isfinite(lon)
    if ok.all():
        return casts
    if ok.sum() < 2:
        raise ValueError("Fewer than 2 casts have any LAT/LON: fill positions from the SADCP first.")
    from scipy.interpolate import interp1d
    f_lat = interp1d(t[ok], lat[ok], fill_value="extrapolate")
    f_lon = interp1d(t[ok], lon[ok], fill_value="extrapolate")
    out = []
    for c, good, ti in zip(casts, ok, t):
        if not good:
            c = c.copy()
            c[COLS["lat"]] = float(f_lat(ti))
            c[COLS["lon"]] = float(f_lon(ti))
        out.append(c)
    print(f"{(~ok).sum()} cast(s) had no LAT/LON -> position interpolated in time from neighbouring casts")
    return out
 
# ---------------------------------------------------------------------------
# TEOS-10 + pressure binning -> 2-D section (pressure x cast)
# ---------------------------------------------------------------------------
def fill_vertical_gaps(a, max_gap):
    """Linearly fill NaN gaps of up to `max_gap` bins inside each column (no extrapolation)."""
    if not max_gap:
        return a
    out = a.copy()
    for j in range(a.shape[1]):
        s = pd.Series(a[:, j])
        out[:, j] = s.interpolate(limit=max_gap, limit_area="inside").values
    return out
 
 
def build_section(casts, dp=1.0, pmax=None, max_gap_dbar=10):
    pmax = pmax or max(c[COLS["p"]].max() for c in casts)
    p_edges = np.arange(0, pmax + dp, dp)
    p_grid = 0.5 * (p_edges[1:] + p_edges[:-1])
 
    n, m = len(p_grid), len(casts)
    SA, CT, sig0 = (np.full((n, m), np.nan) for _ in range(3))
    lon, lat, time = np.zeros(m), np.zeros(m), []
 
    for j, c in enumerate(casts):
        lon[j], lat[j] = c[COLS["lon"]].mean(), c[COLS["lat"]].mean()
        time.append(pd.to_datetime(c[COLS["time"]]).iloc[0])
        p = c[COLS["p"]].values
        sa = gsw.SA_from_SP(c[COLS["sp"]].values, p, lon[j], lat[j])
        ct = gsw.CT_from_t(sa, c[COLS["t"]].values, p)
        idx = np.digitize(p, p_edges) - 1
        ok = (idx >= 0) & (idx < n)
        for arr, val in [(SA, sa), (CT, ct)]:
            s = pd.Series(val[ok]).groupby(idx[ok]).mean()
            arr[s.index.values, j] = s.values
    gap = int(round(max_gap_dbar / dp))
    SA, CT = fill_vertical_gaps(SA, gap), fill_vertical_gaps(CT, gap)
    sig0 = gsw.sigma0(SA, CT)
    return dict(p=p_grid, SA=SA, CT=CT, sig0=sig0, lon=lon, lat=lat, time=pd.DatetimeIndex(time))
 
 
def along_track_x(lon, lat, lon0=None, lat0=None, direction=None):
    """
    Project positions onto the straight line from first to last cast (metres).
    Returns x and the unit vector (east, north) of +x. Using a straight-line
    projection keeps x monotonic and is what the 2-D PV formula assumes.
    """
    lon0 = lon[0] if lon0 is None else lon0
    lat0 = lat[0] if lat0 is None else lat0
    R = 6371e3
    ex = np.deg2rad(lon - lon0) * R * np.cos(np.deg2rad(lat0))
    ny = np.deg2rad(lat - lat0) * R
    if direction is None:
        direction = np.array([ex[-1], ny[-1]]) / np.hypot(ex[-1], ny[-1])
    return ex * direction[0] + ny * direction[1], direction
 
###### SATELLITE FUNCTIONS ######
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
def date_from_name_cmems(f):
    m = re.search(r"(\d{4})-?(\d{2})-?(\d{2})", Path(f).name)
    return pd.Timestamp(f"{m.group(1)}-{m.group(2)}-{m.group(3)}")

def load_adt_cmems(dir_ssh, t0, t1, bbox=None, pad_days=1):
    files = sorted(os.listdir(DIR_SSH))

    # Select cruise dates
    days = pd.date_range(pd.Timestamp(t0) - pd.Timedelta(days=1),
                            pd.Timestamp(t1) + pd.Timedelta(days=1), freq="D")
    tags = {d.strftime("%Y-%m-%d") for d in days} | {d.strftime("%Y%m%d") for d in days}
    picked = [f for f in files if any(t in Path(f).name for t in tags)] or files
    dsets = []
    for f in picked:
        d = xr.open_dataset(dir_ssh / f)

        # fix the time for the daily files 
        if d.sizes.get("time", 1) == 1:
            if "time" in d.dims:
                d = d.isel(time=0, drop=True)
            elif "time" in d.coords:
                d = d.drop_vars("time")
            d = d.expand_dims(time=[date_from_name_cmems(f)])

        # load and close 
        dsets.append(d.load())
        d.close()

    # combine into 1 dataset along time axis
    ds = xr.concat(dsets, "time")
    ds = ds.assign_coords(time=np.asarray(ds.time.values, dtype="datetime64[ns]")).sortby("time")
    ds = ds.isel(time=~pd.Index(ds.time.values).duplicated())
    ds = ds.assign_coords(longitude=xr.where(ds.longitude > 180, ds.longitude - 360, ds.longitude)).sortby(["latitude", "longitude"])
    if bbox is not None:
        x0, x1, y0, y1 = bbox
        ds = ds.sel(longitude=slice(x0, x1), latitude=slice(y0, y1))
    return ds.sel(time=slice(days[0], days[-1] + pd.Timedelta(hours=23)))
ds_cmems = load_adt_cmems(DIR_SSH, CRUISE_START, CRUISE_END,bbox=extent)


###### HYDROGRAPHY FUNCTIONS ######
def determine_sigma(S, T, press, lon, lat):
    # Convert practical salinity to absolute salinity
    SA = gsw.SA_from_SP(S, press, lon, lat)  # Longitude and Latitude needed

    # Convert in situ temperature to conservative temperature
    CT = gsw.CT_from_t(SA, T, press)

    # Calculate in situ density (kg/m³)
    rho = gsw.rho(SA, CT, press)

    # Calculate potential density anomaly referenced to 0 dbar
    sigma0 = gsw.sigma0(SA, CT)  # Potential density anomaly (kg/m³ - 1000)

    return sigma0

###### PLOTTING FUNCTIONS ######
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
 
def plot_section(sec, title="", pmax=None, sig_levels=None, flip_x="auto", savepath=None):
    xkm = sec["x"] / 1e3
    p = sec["p"]
    sig_levels = sig_levels if sig_levels is not None else np.arange(
        np.floor(np.nanmin(sec["sig0"]) * 4) / 4, np.nanmax(sec["sig0"]), 0.25)
 
    fig, axs = plt.subplots(4, 1, figsize=(9, 12), sharex=True, sharey=True,
                            constrained_layout=True)
 
    def panel(ax, data, cmap, label, vmin=None, vmax=None):
        vmin = np.nanpercentile(data, 1) if vmin is None else vmin
        vmax = np.nanpercentile(data, 99) if vmax is None else vmax
        pc = ax.pcolormesh(xkm, p, data, cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
        cs = ax.contour(xkm, p, sec["sig0"], levels=sig_levels, colors="k", linewidths=0.5)
        ax.clabel(cs, cs.levels[::2], fontsize=7, fmt="%.2f")
        fig.colorbar(pc, ax=ax, label=label, pad=0.01)
        ax.set_ylabel("Pressure (dbar)")
        return pc
 
    panel(axs[0], sec["CT"], CMAPS["CT"], r"$\Theta$ (°C)")
    axs[0].set_title("Conservative Temperature")
    panel(axs[1], sec["SA"], CMAPS["SA"], r"$S_A$ (g kg$^{-1}$)")
    axs[1].set_title("Absolute Salinity")
    panel(axs[2], sec["sig0"], CMAPS["sig"], r"$\sigma_0$ (kg m$^{-3}$)")
    axs[2].set_title(r"Potential density anomaly $\sigma_0$")
 
    # PV normalised by f * N2_ref is hard to read across depth -> show q directly, symmetric
    qplot = sec["q"] * 1e9
    qmax = np.nanpercentile(np.abs(qplot), 95)
    panel(axs[3], qplot, CMAPS["pv"], r"$q$ ($10^{-9}$ s$^{-3}$)", -qmax, qmax)
    # mark where f*q < 0 (symmetric/inertial instability possible)
    unstable = (sec["f"] * sec["q"] < 0).astype(float)
    axs[3].contour(xkm, p, unstable, levels=[0.5], colors="magenta", linewidths=1)
    zeta_txt = "with ADCP ζ" if np.any(sec["zeta"]) else "ζ = 0 (no ADCP)"
    axs[3].set_title(rf"Ertel PV  $q=(f+\zeta)N^2-b_x^2/f$  [{zeta_txt}];  magenta: $fq<0$")
 
    axs[0].invert_yaxis()
    if pmax:
        axs[0].set_ylim(pmax, 0)
    if flip_x == "auto":
        flip_x = auto_flip(sec["xdir"])
    if flip_x:
        axs[0].invert_xaxis()
    axs[-1].set_xlabel("Along-track distance (km)")
    fig.suptitle(title)
    if savepath:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")
    return fig, axs
 
def draw_background(ax, crs, ssh, eddy_files, transect,eddy_source, adt_ticks, date, levels):
    """ADT + eddies + light bathymetry. The transect is drawn on top by the module."""
    c = ax.contourf(ssh["lon"], ssh["lat"], ssh["adt"][0], levels=levels,
                    cmap="RdBu_r", extend="both", alpha=0.85, transform=crs)
    ax.contour(lon_subset, lat_subset, depth_subset, levels=[-4000, -3000, -2000, -1000],
               colors="0.35", linewidths=0.4, alpha=0.5, transform=crs)
    if eddy_files:
        anti, cyc = load_eddies(eddy_files)
        plot_eddies(ax, anti, "darkred", crs)
        plot_eddies(ax, cyc, "navy", crs)
    handles = [mlines.Line2D([], [], color="darkred", label="Anticyclone (max speed)"),
            #    mlines.Line2D([], [], color="darkred", ls="--", label="Anticyclone (outer)"),
               mlines.Line2D([], [], color="navy", label="Cyclone (max speed)"),
            #    mlines.Line2D([], [], color="navy", ls="--", label="Cyclone (outer)"),
               mlines.Line2D([], [], color="k", lw=2.5, label=f"MVP {transect}")]
    leg = ax.legend(handles=handles, fontsize=8, loc="lower left", framealpha=0.9)
    leg.set_zorder(200)
    if eddy_source:
        ax.text(0.99, 0.01, eddy_source, transform=ax.transAxes, fontsize=6, ha="right",
                va="bottom", zorder=200,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
    cb = ax.figure.colorbar(c, ax=ax, orientation="horizontal", shrink=0.8, pad=0.05,
                            ticks=adt_ticks)
    cb.set_label(f"ADT {date} [m]")
 
def mark_section_clean(ax, x_oc, x_ms, x_c, x_nen, z, col, vertical_nen, tag="", label_center=True):
    """Markers for ONE eddy on one section panel: max speed solid, outer dashed (thin),
    TOEddies center black dashed, Nencioli green dashed. NaN centres are skipped.
    tag : short eddy name for the labels, e.g. "C" or "AC"."""
    tag = f"{tag} " if tag else ""
    xt = ax.get_xaxis_transform()
    for xv in x_oc:
        ax.axvline(xv, color=col, ls="--", lw=1.0, alpha=0.8)
    for xv in x_ms:
        ax.axvline(xv, color=col, ls="-", lw=1.6)
    if np.isfinite(x_c):
        ax.axvline(x_c, color="k", ls="--", lw=1.2)
        if label_center:
            ax.text(x_c, 1.01, f"{tag}TOEddies\ncenter ", transform=xt,      # left of the line
                    ha="right", va="bottom", fontsize=6.5)
    if vertical_nen:
        if np.isfinite(x_nen):
            ax.axvline(x_nen, color=NEN_COLOR, ls="--", lw=1.4)
            if label_center:
                ax.text(x_nen, 1.01, f" {tag}Nencioli\n center", color=NEN_COLOR, transform=xt,
                        ha="left", va="bottom", fontsize=6.5)                 # right of the line
    else:
        x_nen = np.asarray(x_nen, float)
        good = np.flatnonzero(np.isfinite(x_nen))
        if good.size:
            ax.plot(x_nen, z, color=NEN_COLOR, ls="--", lw=1.4)
            if label_center:                    # label at the shallowest level in frame
                ax.text(x_nen[good[0]], 1.01, f" {tag}Nencioli\n center", color=NEN_COLOR,
                        transform=xt, ha="left", va="bottom", fontsize=6.5)
 
 
def mark_eddies_on_sections(axs, sec, z, a_lon, a_lat, u, v, edds, min_fixes=5):
    """Draw the markers of every eddy in `edds` (list from eddies_on_transect) on the
    3 MVP panels axs[:3] and the SADCP panel axs[3]. Only what falls inside the section
    (0 .. end km along the line) is drawn.
 
    Nencioli is run separately for each eddy, using ONLY the SADCP fixes inside that
    eddy's outer contour and a search box no larger than the eddy itself, so the
    neighbouring eddy's velocities cannot pull the centre away. A centre found on the
    edge of the search box (= no real minimum) is dropped."""
    xkm, p = sec["x"] / 1e3, sec["p"]
    to_sn = make_to_sn(sec)
    x_lo, x_hi = np.nanmin(xkm), np.nanmax(xkm)
    a_lon, a_lat = np.asarray(a_lon, float), np.asarray(a_lat, float)
 
    def in_frame(x):
        x = np.atleast_1d(np.asarray(x, float))
        return np.isfinite(x) & (x >= x_lo) & (x <= x_hi)
 
    for e in edds:
        tag = "C" if e.get("cyc", True) else "AC"
        col = "tab:blue" if e.get("cyc", True) else "tab:red"
        clon, clat = e["center"]
 
        x_oc = crossings_on_line(e["outer"], to_sn)
        x_ms = crossings_on_line(e["maxspeed"], to_sn)
        x_oc, x_ms = x_oc[in_frame(x_oc)], x_ms[in_frame(x_ms)]
 
        x_c = float(to_sn(clon, clat)[0])
        if not in_frame(x_c)[0]:
            print(f"  {tag}: TOEddies center at {x_c:.0f} km, outside section "
                  f"({x_lo:.0f}-{x_hi:.0f} km) -> not marked")
            x_c = np.nan
 
        # --- Nencioli, this eddy only (layer centre + one centre per depth) ---
        x_nen, x_nen_z = np.nan, np.full(len(z), np.nan)
        nen = nencioli_eddy(e, dict(lon=a_lon, lat=a_lat, z=z, u=u, v=v),
                            min_fixes=min_fixes, by_depth=True)
        if nen is None:
            print(f"  {tag}: too few SADCP fixes inside the outer contour -> no Nencioli")
        else:
            xz = to_sn(nen["lon_z"], nen["lat_z"])[0]
            x_nen_z = np.where(in_frame(xz) & ~nen["edge_z"], xz, np.nan)
            xl = float(to_sn(nen["lon"], nen["lat"])[0])
            if nen["edge"]:
                print(f"  {tag}: Nencioli minimum on the edge of the search box -> no centre")
            elif not in_frame(xl)[0]:
                print(f"  {tag}: Nencioli center at {xl:.0f} km, outside section -> not marked")
            else:
                x_nen = xl

        for ax in axs[:3]:                                   # MVP panels: vertical lines
            mark_section_clean(ax, x_oc, x_ms, x_c, x_nen, p, col, vertical_nen=True, tag=tag)
        mark_section_clean(axs[3], x_oc, x_ms, x_c, x_nen_z, z, col, vertical_nen=False, tag=tag)
 
 
def plot_section_with_map(sec, adcp, map_draw=None, map_args=None, title="", pmax=None, sig_levels=None,
                          map_extent=None, map_pad=(2.5, 2.0), vel_lim=None, edd=None, savepath=None,
                          flip_x="auto", woa_lon=None, woa_lat=None):
    """
    Layout: 2 x 2 sections on the left
                (a) Theta        (b) S_A
                (c) sigma0       (d) SADCP cross-track velocity
            small map on the right.
 
    map_draw   : function(ax, crs, *map_args) that draws the background (ADT, eddies, bathymetry),
                 e.g. draw_background with map_args = [ssh, eddy_files, transect, eddy_source,
                 adt_ticks, date, levels].
    vel_lim    : colour limit for velocity (symmetric). Default = 98th pct of |v|.
    pmax       : max pressure/depth shown on all 4 panels. Default = deepest data of each panel.
    map_extent : [lon_min, lon_max, lat_min, lat_max]; default = transect +/- map_pad (deg).
    edd        : list of eddies from eddies_on_transect (or a single dict); all are marked.
    flip_x     : "auto" -> the x axis is oriented like the map (W left / E right for a zonal
                 transect, S left / N right for a meridional one), so the ship moves across
                 the panels the same way it moves on the map. True -> 0 km on the right,
                 False -> 0 km on the left.
    """
    import matplotlib.patheffects as pe
    from matplotlib.gridspec import GridSpec
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        ccrs = None
 
    lon, lat = sec["lon"], sec["lat"]
    xkm, p = sec["x"] / 1e3, sec["p"]
    if map_extent is None:
        map_extent = [np.nanmin(lon) - map_pad[0], np.nanmax(lon) + map_pad[0],
                      np.nanmin(lat) - map_pad[1], np.nanmax(lat) + map_pad[1]]
    sig_levels = sig_levels if sig_levels is not None else np.arange(
        np.floor(np.nanmin(sec["sig0"]) * 4) / 4, np.nanmax(sec["sig0"]), 0.25)
 
    fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    # left = 2 x 2 sections, right = small map. Separate sub-grids so the map
    # (and its colourbar) does not change the size of the section panels.
    outer = GridSpec(1, 2, figure=fig, width_ratios=[4, 1])
    gs_sec = outer[0].subgridspec(2, 2)
    gs_map = outer[1].subgridspec(2, 1, height_ratios=[1, 0.6])   # map in the top part
 
    # ======================= sections (left, 2 x 2) =======================
    # x is shared by all 4 panels; y is NOT shared, so each panel gets its own
    # (inverted) depth range: MVP pressure for (a)-(c), SADCP depth for (d).
    ax0 = fig.add_subplot(gs_sec[0, 0])
    axs = [ax0,
           fig.add_subplot(gs_sec[0, 1], sharex=ax0),
           fig.add_subplot(gs_sec[1, 0], sharex=ax0),
           fig.add_subplot(gs_sec[1, 1], sharex=ax0)]
 
    specs = [(sec["CT"], CMAPS["CT"], r"$\Theta$ (°C)", "(a) Conservative Temperature"),
             (sec["SA"], CMAPS["SA"], r"$S_A$ (g kg$^{-1}$)", "(b) Absolute Salinity"),
             (sec["sig0"], CMAPS["sig"], r"$\sigma_0$ (kg m$^{-3}$)", r"(c) Potential density $\sigma_0$")]
    for ax, (data, cmap, lab, ttl) in zip(axs[:3], specs):
        pc = ax.pcolormesh(xkm, p, data, cmap=cmap, shading="nearest",
                           vmin=np.nanpercentile(data, 1), vmax=np.nanpercentile(data, 99))
        fig.colorbar(pc, ax=ax, label=lab, pad=0.01)
        ax.set_title(ttl, loc="left", fontsize=11, pad=14)
 
    # --- (d) SADCP cross-track velocity ---
    axc = axs[3]
    a_lon, a_lat, z, u, v = adcp["lon"], adcp["lat"], adcp["z"], adcp["u"], adcp["v"]
 
    # same origin (first cast) and same direction (first -> last cast) as sec["x"]
    s_m, _ = along_track_x(np.asarray(a_lon), np.asarray(a_lat),
                           sec["lon"][0], sec["lat"][0], sec["xdir"])
    s = s_m / 1e3
    ex = sec["xdir"]
    # only the SADCP profiles along the MVP section (matters when the section is cut)
    x_lo, x_hi = np.nanmin(xkm), np.nanmax(xkm)
    on = np.isfinite(s) & (s >= x_lo - DS_KM) & (s <= x_hi + DS_KM)
    a_lon, a_lat, u, v, s = (np.asarray(a)[on] for a in (a_lon, a_lat, u, v, s))
    vn = -u * ex[1] + v * ex[0]        # positive = to the LEFT of the 0 km -> end direction
 
    edges = np.arange(np.floor(np.nanmin(s)), np.ceil(np.nanmax(s)) + DS_KM, DS_KM)
    xg, V, _ = bin_columns(s, vn, edges)
    V = fill_z_gaps(V, z)
 
    if vel_lim is None:
        vals = np.abs(V.ravel())
        vlim = max(0.1, np.round(np.nanpercentile(vals[np.isfinite(vals)], 98) * 10) / 10)
    else:
        vlim = vel_lim
    fill_lev = np.linspace(-vlim, vlim, N_FILL + 1)
    cf = axc.contourf(xg, z, V.T, levels=fill_lev, cmap=cmocean.cm.balance, extend="both")
    nstep = int(np.ceil(vlim / C_LINE_STEP)) + 1
    line_lev = np.round(np.arange(-nstep, nstep + 1) * C_LINE_STEP, 3)
    cs = axc.contour(xg, z, V.T, levels=line_lev, colors="k", linewidths=0.6)
    axc.clabel(cs, cs.levels, fmt="%.1f", fontsize=6, inline=True, inline_spacing=2)
    cbc = fig.colorbar(cf, ax=axc, pad=0.01)
    cbc.set_ticks(fill_lev[1::3])
    cbc.ax.set_yticklabels([f"{0.0 if abs(t) < 1e-9 else t:.2f}" for t in fill_lev[1::3]])
    cbc.set_label(r"v$_\perp$ (m s$^{-1}$)")
    axc.set_title(fr"(d) SADCP: v$_\perp$ {DS_KM:.0f} km bins", loc="left", fontsize=11, pad=14)
 
    # ---- TOEddies / Nencioli markers on all 4 panels (every eddy crossed) ----
    # edd: one dict (old eddy_on_transect) or a list of dicts (eddies_on_transect)
    edds = [edd] if isinstance(edd, dict) else list(edd or [])
    if edds:
        mark_eddies_on_sections(axs, sec, z, a_lon, a_lat, u, v, edds)
 
    # ---- common decorations on all 4 sections ----
    for ax in axs:
        cs = ax.contour(xkm, p, sec["sig0"], levels=sig_levels, colors="k", linewidths=0.5)
        ax.clabel(cs, cs.levels[::2], fontsize=7, fmt="%.2f")
        ax.plot(xkm, np.zeros_like(xkm), "|", color="k", ms=6, clip_on=False)  # cast positions
 
    # y axes: each panel inverted on its own (surface at the top)
    p_bot = pmax if pmax else np.nanmax(p)
    z_bot = pmax if pmax else np.nanmax(z)
    for ax in axs[:3]:
        ax.set_ylim(p_bot, 0)
    axc.set_ylim(z_bot, 0)
    for ax in axs[:3]:
        ax.set_ylabel("Pressure (dbar)")
    axc.set_ylabel("Depth (m)")
 
    # x axis: shared, labels only on the bottom row
    axs[0].set_xlim(np.nanmin(xkm), np.nanmax(xkm))
    if flip_x == "auto":
        flip_x = auto_flip(sec["xdir"])
    if flip_x:
        axs[0].invert_xaxis()                     # applies to all 4 (sharex)
    for ax in axs[:2]:                            # top row: ticks, no labels
        ax.tick_params(labelbottom=False)
    for ax in axs[2:]:                            # bottom row: tick labels + axis label
        ax.tick_params(labelbottom=True)
        ax.set_xlabel("Along-track distance from start (km)")
 
    # ======================= map (small, right) =======================
    if ccrs is not None:
        crs = ccrs.PlateCarree()
        proj = ccrs.Mercator(central_longitude=np.mean(map_extent[:2]))
        axm = fig.add_subplot(gs_map[0], projection=proj)
        axm.set_extent(map_extent, crs=crs)
        tr = dict(transform=crs)
    else:                                   # fallback without cartopy
        crs, tr = None, {}
        axm = fig.add_subplot(gs_map[0])
        axm.set_xlim(map_extent[:2]); axm.set_ylim(map_extent[2:])
        axm.set_aspect(1 / np.cos(np.deg2rad(np.mean(map_extent[2:]))))
 
    if map_draw is not None:
        map_draw(axm, crs, *(map_args or []))
 
    if ccrs is not None:
        axm.add_feature(cfeature.LAND, facecolor="0.8", zorder=50)
        axm.coastlines(resolution="50m", linewidth=0.6, zorder=51)
        gl = axm.gridlines(crs=crs, draw_labels=True, linewidth=0.4, color="0.5",
                           alpha=0.4, linestyle="--")
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {"size": 8}
    else:
        axm.grid(lw=0.4, alpha=0.4, ls="--")
        axm.tick_params(labelsize=8)
 
    halo = [pe.Stroke(linewidth=4.5, foreground="white"), pe.Normal()]
    axm.plot(lon, lat, "-", color="k", lw=2, path_effects=halo, zorder=100, **tr)
    if woa_lon is not None and len(woa_lon):
        axm.scatter(woa_lon, woa_lat, s=6, color="orange", edgecolor="k", linewidth=0.2,
                    zorder=101, label=f"WOA grid points ({len(woa_lon)})", **tr)
    h, l = axm.get_legend_handles_labels()
    leg = axm.get_legend()
    if leg is not None:
        axm.legend(handles=leg.legend_handles + h, fontsize=8, loc="lower left",
                       framealpha=0.9).set_zorder(200)
    for i, lab in [(0, "0 km"), (-1, f"{xkm[-1]:.0f} km")]:
        axm.plot(lon[i], lat[i], marker="o" if i == 0 else "s", ms=6, mfc="yellow",
                 mec="k", mew=1, zorder=102, **tr)
        axm.annotate(lab, (lon[i], lat[i]), xytext=(6, 0), textcoords="offset points",
                     fontsize=8, fontweight="bold", va="center", zorder=103,
                     path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
                     xycoords=crs._as_mpl_transform(axm) if crs is not None else "data")
    t0 = sec["time"][0].strftime("%Y-%m-%d %H:%M")
    t1 = sec["time"][-1].strftime("%Y-%m-%d %H:%M")
    axm.set_title(f"{t0}\nto {t1} UTC", fontsize=9)
 
    fig.suptitle(title, fontsize=14, fontweight="bold")
    if savepath:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")
    return fig, [axm] + axs

###### EDDY CHARACTERISTICS (WOA anomalies, zeta/f, EPV core, eddy metrics) ######
# Same methods and parameters as anom_woa.ipynb and epv_calc.ipynb, as functions.
from scipy.ndimage import uniform_filter1d
from scipy import ndimage

DIR_WOA = BASE / "data/woa"
WOA_RADIUS_KM = 100       # WOA grid points within this distance of the eddy centre
RHO0 = 1025.0             # kg m^-3
DZ_EPV = 8.0              # m, depth bins of the EPV grid
ZMAX_EPV = 1000.0         # m, bottom of the EPV grid
SMOOTH_X_V = 3            # running mean on v_perp along x (bins)
SMOOTH_Z_V = 3            # running mean on v_perp in depth (SADCP levels)
SMOOTH_Z_B = 3            # running mean on b in depth (DZ_EPV bins)
MAX_XGAP_KM = 10          # fill MVP gaps along x up to this size (km)
REF_ZMIN = 10.0           # m: background levels shallower than this are ignored
MLD_DSIG = 0.03           # kg m^-3: background mixed-layer criterion
WOA_N2_SMOOTH_M = 150     # m, running mean on N^2_ref for WOA
EPV_X_MIN = 1e-13         # s^-3: R masked where |EPV_x| is below this
BETA = 40                 # eddy core: |dEPV_z| > BETA |EPV_x|
CORE_SMOOTH = 3           # bins, smoothing before the core test
DEPV_MIN = 1e-11          # s^-3: |dEPV_z| must also exceed this
CORE_ZMIN = 60            # m: automatic seed searched below this depth
CORE_SIGN = True          # require dEPV_z of the eddy's sign


def nan_smooth(a, size, axis):
    """Running mean that ignores NaNs (NaNs stay NaN)."""
    if size <= 1:
        return a
    w = np.isfinite(a).astype(float)
    num = uniform_filter1d(np.where(w > 0, a, 0.0), size, axis=axis, mode="nearest")
    den = uniform_filter1d(w, size, axis=axis, mode="nearest")
    out = num / np.where(den > 0, den, np.nan)
    out[w == 0] = np.nan
    return out


# ---------------- WOA reference and anomalies on isopycnals ----------------
def woa_reference_around(lon_c, lat_c, radius_km=WOA_RADIUS_KM, dir_woa=DIR_WOA):
    """Mean WOA (July) profile of the grid points within radius_km of (lon_c, lat_c).
    -> ref = dict(p, CT, SA, sig0) on the WOA levels with data, and the lon/lat of the points used."""
    f_t = dir_woa / "temp" / sorted(os.listdir(dir_woa / "temp"))[0]
    f_s = dir_woa / "salt" / sorted(os.listdir(dir_woa / "salt"))[0]
    with nc.Dataset(f_t) as dt, nc.Dataset(f_s) as dsl:
        lon2d, lat2d = np.meshgrid(np.array(dt["lon"][:]), np.array(dt["lat"][:]))
        temp = np.ma.filled(dt["t_an"][:], np.nan)
        sal = np.ma.filled(dsl["s_an"][:], np.nan)
        depth = np.array(dt["depth"][:], float)
    dist = np.hypot((lat2d - lat_c) * 111.0, (lon2d - lon_c) * 111.0 * np.cos(np.deg2rad(lat_c)))
    mask = (dist <= radius_km) & np.isfinite(temp[0, 0]) & np.isfinite(sal[0, 0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        t = np.squeeze(np.nanmean(temp[..., mask], axis=-1))
        s = np.squeeze(np.nanmean(sal[..., mask], axis=-1))
    lo, la = lon2d[mask].mean(), lat2d[mask].mean()
    p = gsw.p_from_z(-depth, la)
    SA = gsw.SA_from_SP(s, p, lo, la)
    CT = gsw.CT_from_t(SA, t, p)
    ok = np.isfinite(SA) & np.isfinite(CT)
    ref = dict(p=p[ok], CT=CT[ok], SA=SA[ok], sig0=gsw.sigma0(SA, CT)[ok], file=f_t.name)
    return ref, lon2d[mask], lat2d[mask]


def interp_on_rho(sigma, var, ind_rho):
    """var(sigma) -> var on the density levels ind_rho (NaN outside the sigma range)."""
    ok = np.isfinite(sigma) & np.isfinite(var)
    if ok.sum() < 2:
        return np.full(len(ind_rho), np.nan)
    order = np.argsort(sigma[ok])
    s, idx = np.unique(sigma[ok][order], return_index=True)
    if s.size < 2:
        return np.full(len(ind_rho), np.nan)
    return np.interp(ind_rho, s, var[ok][order][idx], left=np.nan, right=np.nan)


def iso_anomaly(ref, sec, var="CT", pas=0.01):
    """Anomaly of sec[var] (CT or SA) relative to ref ON ISOPYCNALS (Barabinot et al. 2025),
    put back at the pressure of each isopycnal in each cast -> same shape as sec["CT"].
    var="sig0": density anomaly at constant pressure instead."""
    p, sig = sec["p"], sec["sig0"]
    if var == "sig0":
        return sig - np.interp(p, ref["p"], ref["sig0"], left=np.nan, right=np.nan)[:, None]
    ind_rho = np.arange(22, 28, pas)
    ref_i = interp_on_rho(ref["sig0"], ref[var], ind_rho)
    out = np.full(sig.shape, np.nan)
    for i in range(sig.shape[1]):
        a_rho = interp_on_rho(sig[:, i], sec[var][:, i], ind_rho) - ref_i
        p_rho = interp_on_rho(sig[:, i], p, ind_rho)
        good = np.isfinite(a_rho) & np.isfinite(p_rho)
        if good.sum() < 2:
            continue
        o = np.argsort(p_rho[good])
        d, idx = np.unique(p_rho[good][o], return_index=True)
        if d.size >= 2:
            out[:, i] = np.interp(p, d, a_rho[good][o][idx], left=np.nan, right=np.nan)
    return out


# ---------------- zeta/f and Ertel PV on the section ----------------
def section_vorticity(sec, adcp, dx_km=DS_KM, smooth_x=SMOOTH_X_V, smooth_z=SMOOTH_Z_V):
    """SADCP cross-track velocity (+ to the LEFT of +x) in dx_km bins on the section axis,
    zeta = dv_perp/dx, f and zeta/f. Only the SADCP profiles along the (possibly cut) section."""
    z = np.asarray(adcp["z"], float)
    lon, lat = np.asarray(adcp["lon"], float), np.asarray(adcp["lat"], float)
    s = along_track_x(lon, lat, sec["lon"][0], sec["lat"][0], sec["xdir"])[0] / 1e3
    xkm = sec["x"] / 1e3
    on = np.isfinite(s) & (s >= np.nanmin(xkm) - dx_km) & (s <= np.nanmax(xkm) + dx_km)
    ex = sec["xdir"]
    vp = -np.asarray(adcp["u"])[on] * ex[1] + np.asarray(adcp["v"])[on] * ex[0]
    edges = np.arange(np.floor(np.nanmin(s[on])), np.ceil(np.nanmax(s[on])) + dx_km, dx_km)
    xg, V, cnt = bin_columns(s[on], vp, edges)
    _, LATb, _ = bin_columns(s[on], lat[on][:, None], edges)
    V = fill_z_gaps(V, z)
    V = nan_smooth(nan_smooth(V, smooth_x, axis=0), smooth_z, axis=1)
    zeta = np.gradient(V, xg * 1e3, axis=0)
    lat_b = LATb[:, 0]
    lat_b = np.where(np.isfinite(lat_b), lat_b, np.nanmean(lat_b))
    f = 2 * OMEGA * np.sin(np.deg2rad(lat_b))
    return dict(x=xg, z=z, v_perp=V, zeta=zeta, f=f, ro=zeta / f[:, None], n_fixes=cnt, dx_km=dx_km)


def _interp_rows(A, z_in, z_out):
    out = np.full((A.shape[0], len(z_out)), np.nan)
    for i, row in enumerate(A):
        ok = np.isfinite(row)
        if ok.sum() >= 2:
            out[i] = np.interp(z_out, z_in[ok], row[ok], left=np.nan, right=np.nan)
    return out


def section_epv(sec, vort, dz=DZ_EPV, zmax=ZMAX_EPV):
    """Ertel PV on the vorticity x bins and dz depth bins (epv_calc.ipynb):
    EPV_z = (f + zeta) N^2, EPV_x = -(dv_perp/dz)(db/dx), b from the MVP sigma0."""
    xg, z_a, dx = vort["x"], vort["z"], vort["dx_km"]
    edges_x = np.r_[xg - dx / 2, xg[-1] + dx / 2]
    lat0 = np.nanmean(sec["lat"])
    depth_mvp = -gsw.z_from_p(sec["p"], lat0)
    has_v = np.isfinite(vort["v_perp"]).any(axis=0)
    z_top = np.nanmin(z_a[has_v])
    z_bot = min(zmax, np.nanmax(z_a[has_v]), np.nanmax(depth_mvp[np.isfinite(sec["sig0"]).any(axis=1)]))
    z_edges = np.arange(z_top - dz / 2, z_bot + dz / 2, dz)
    z = 0.5 * (z_edges[1:] + z_edges[:-1])

    V = _interp_rows(vort["v_perp"], z_a, z)
    zeta = np.gradient(V, xg * 1e3, axis=0)
    dvdz = -np.gradient(V, z, axis=1)                       # z up
    b_sec = -G * (sec["sig0"] + 1000.0 - RHO0) / RHO0
    _, Bx, _ = bin_columns(sec["x"] / 1e3, b_sec.T, edges_x)
    _, Bxz, _ = bin_columns(depth_mvp, Bx.T, z_edges)
    B = pd.DataFrame(Bxz.T).interpolate(axis=0, limit=int(MAX_XGAP_KM // dx),
                                        limit_area="inside").values
    B = nan_smooth(nan_smooth(B, SMOOTH_X_V, axis=0), SMOOTH_Z_B, axis=1)
    N2 = -np.gradient(B, z, axis=1)
    dbdx = np.gradient(B, xg * 1e3, axis=0)
    f = vort["f"]
    EPV_z = (f[:, None] + zeta) * N2
    EPV_x = -dvdz * dbdx
    return dict(x=xg, z=z, f=f, v_perp=V, zeta=zeta, dvdz=dvdz, b=B, N2=N2, dbdx=dbdx,
                EPV=EPV_z + EPV_x, EPV_z=EPV_z, EPV_x=EPV_x, lat0=lat0,
                sig0=-(B * RHO0 / G) + RHO0 - 1000.0)


def epv_background(ref, lat0, dz=DZ_EPV):
    """Background b, sigma0, N^2 on its own depth grid (CTD: dense; WOA: coarse levels)."""
    d = -gsw.z_from_p(ref["p"], lat0)
    b = -G * (ref["sig0"] + 1000.0 - RHO0) / RHO0
    o = np.argsort(d)
    d, b = d[o], b[o]
    keep = d >= REF_ZMIN
    d, b = d[keep], b[keep]
    sig_raw = -(b * RHO0 / G) + RHO0 - 1000.0
    below = np.where(sig_raw >= sig_raw[0] + MLD_DSIG)[0]
    mld = d[below[0]] if below.size else d[-1]
    N2_ml = -(np.interp(mld, d, b) - b[0]) / (mld - d[0])
    edges = np.arange(REF_ZMIN, d.max() + dz, dz)
    zr = 0.5 * (edges[1:] + edges[:-1])
    if np.nanmedian(np.diff(d)) <= dz:                       # dense (CTD)
        _, bb, _ = bin_columns(d, b[:, None], edges)
        bb = nan_smooth(bb[:, 0][None, :], SMOOTH_Z_B, axis=1)[0]
        N2 = nan_smooth(-np.gradient(bb, zr)[None, :], SMOOTH_Z_B, axis=1)[0]
    else:                                                    # coarse (WOA)
        bb = np.interp(zr, d, b, left=np.nan, right=np.nan)
        N2 = np.interp(zr, 0.5 * (d[1:] + d[:-1]), -np.diff(b) / np.diff(d), left=np.nan, right=np.nan)
        N2 = nan_smooth(N2[None, :], max(1, int(round(WOA_N2_SMOOTH_M / dz))), axis=1)[0]
    return dict(z=zr, b=bb, sig0=-(bb * RHO0 / G) + RHO0 - 1000.0, N2=N2, mld=mld, N2_ml=N2_ml)


def epv_anomaly_iso(epv, ref, outcrop="nan"):
    """dEPV_z = EPV_z - f N^2_ref(sigma0) (background N^2 on the same isopycnal), R = dEPV_z/EPV_x."""
    prof = epv_background(ref, epv["lat0"])
    ok = np.isfinite(prof["sig0"]) & np.isfinite(prof["N2"])
    s, i = np.unique(np.maximum.accumulate(prof["sig0"][ok]), return_index=True)
    left = prof["N2_ml"] if outcrop == "surface" else np.nan
    N2r = np.interp(epv["sig0"], s, prof["N2"][ok][i], left=left, right=np.nan)
    f2 = epv["f"][:, None]
    dEPV_z = epv["EPV_z"] - f2 * N2r
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(np.abs(epv["EPV_x"]) > EPV_X_MIN, dEPV_z / epv["EPV_x"], np.nan)
    return dict(prof=prof, N2r=N2r, EPV_ref=f2 * N2r, dEPV_z=dEPV_z, R=R, sig_ref_min=s[0])


def epv_core_mask(epv, an, eddy, to_sn, beta=BETA, seed=None):
    """Eddy core (Barabinot et al. sect. 4.5): connected region around the seed where
    dEPV_z (of the eddy's sign) > beta |EPV_x|. Seed = TOEddies centre column, depth of
    max dEPV_z below CORE_ZMIN. -> (mask (nx, nz), seed (x_km, z_m))."""
    sgn = -1.0 if eddy["cyc"] else 1.0                       # SH: AC -> dEPV_z > 0
    dE_s = nan_smooth(nan_smooth(an["dEPV_z"], CORE_SMOOTH, 0), CORE_SMOOTH, 1)
    dE = sgn * dE_s if CORE_SIGN else np.abs(dE_s)
    Ex = nan_smooth(nan_smooth(np.abs(epv["EPV_x"]), CORE_SMOOTH, 0), CORE_SMOOTH, 1)
    crit = np.isfinite(dE) & np.isfinite(Ex) & (dE > beta * Ex) & (dE > DEPV_MIN)
    if seed is None:
        x_c = float(to_sn(*eddy["center"])[0])
        ix = int(np.argmin(np.abs(epv["x"] - x_c)))
        col = np.where(epv["z"] >= CORE_ZMIN, sgn * an["dEPV_z"][ix], np.nan)
        if not np.isfinite(col).any():
            return np.zeros_like(crit), (epv["x"][ix], np.nan)
        seed = (epv["x"][ix], epv["z"][int(np.nanargmax(col))])
    ix = int(np.argmin(np.abs(epv["x"] - seed[0])))
    iz = int(np.argmin(np.abs(epv["z"] - seed[1])))
    lab, _ = ndimage.label(crit, structure=np.ones((3, 3)))
    if lab[ix, iz] == 0:
        return np.zeros_like(crit), seed
    return ndimage.binary_fill_holes(lab == lab[ix, iz]), seed


def core_on_section(core, epv, sec):
    """EPV core mask (nx, nz on the EPV grid) -> boolean mask on the MVP section grid
    (pressure x cast, like sec["CT"]): True where the cast's x bin and the level's depth
    bin are inside the core."""
    dx = np.median(np.diff(epv["x"]))
    dz = np.median(np.diff(epv["z"]))
    ix = np.rint((sec["x"] / 1e3 - epv["x"][0]) / dx).astype(int)
    depth = -gsw.z_from_p(sec["p"], epv["lat0"])
    iz = np.rint((depth - epv["z"][0]) / dz).astype(int)
    okx = (ix >= 0) & (ix < core.shape[0])
    okz = (iz >= 0) & (iz < core.shape[1])
    m = np.zeros((len(sec["p"]), len(sec["x"])), bool)
    m[np.ix_(okz, okx)] = core[np.ix_(ix[okx], iz[okz])].T
    return m


# ---------------- Nencioli centre of one eddy ----------------
def nencioli_eddy(eddy, adcp, min_fixes=5, layer=NEN_LAYER, by_depth=False):
    """Nencioli centre of ONE eddy from the SADCP fixes inside its outer contour, in a search box
    no larger than the eddy (so a neighbouring eddy cannot pull the centre away).
    adcp: dict(lon, lat, z, u, v). Returns None if < min_fixes fixes inside, else
    dict(lon, lat, edge, mtg) for the `layer` mean, n_fixes, search_km, and with by_depth=True
    also lon_z, lat_z, edge_z (one centre per SADCP level). edge = minimum on the box edge
    (no real minimum)."""
    a_lon, a_lat = np.asarray(adcp["lon"], float), np.asarray(adcp["lat"], float)
    z = np.asarray(adcp["z"], float)
    inside = MplPath(eddy["outer"]).contains_points(np.c_[a_lon, a_lat])
    inside &= np.isfinite(a_lon) & np.isfinite(a_lat)
    if inside.sum() < min_fixes:
        return None
    clon, clat = eddy["center"]
    ox, oy = to_xy(eddy["outer"][:, 0], eddy["outer"][:, 1], clon, clat)
    search = min(NEN_SEARCH_KM, float(np.nanmax(np.hypot(ox, oy))))   # outer radius
    lo, la = a_lon[inside], a_lat[inside]
    uu, vv = np.asarray(adcp["u"])[inside], np.asarray(adcp["v"])[inside]
    lay = (z >= layer[0]) & (z <= layer[1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ul, vl = np.nanmean(uu[:, lay], axis=1), np.nanmean(vv[:, lay], axis=1)
    nl = nencioli_center(lo, la, ul[:, None], vl[:, None], clon, clat, search_km=search)
    out = dict(lon=float(nl["lon"][0]), lat=float(nl["lat"][0]), edge=bool(nl["edge"][0]),
               mtg=float(nl["mtg"][0]), n_fixes=int(inside.sum()), search_km=search)
    if by_depth:
        nd = nencioli_center(lo, la, uu, vv, clon, clat, search_km=search)
        out.update(lon_z=nd["lon"], lat_z=nd["lat"], edge_z=nd["edge"])
    return out


def contour_radius_km(xy, lon0, lat0):
    """Equivalent radius sqrt(area / pi) of a closed lon/lat contour, km."""
    x, y = to_xy(xy[:, 0], xy[:, 1], lon0, lat0)
    area = 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    return float(np.sqrt(area / np.pi))
