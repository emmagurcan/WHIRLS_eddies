"""
Orthogonal-velocity figures of the WHIRLS eddies from SADCP and LADCP.

Two ways of looking at an eddy, chosen per section in EDDIES below:

mode = "section"  (straight crossing, e.g. the cyclone CE18)
    A straight line is fitted through the ship positions of the time window. Every
    measurement is projected on it: x = distance along the line, and the velocity
    component perpendicular to the line (positive to the LEFT of the ship's direction of
    travel) is plotted against distance and depth, like Fig. 3 of Subirade et al. (2023).
    The TOEddies center and the max-speed radius are projected on the line.

mode = "radial"   (ship going around the eddy, e.g. the anticyclone AE11)
    A straight line makes no sense. Every measurement is put in the eddy's moving frame
    (relative to the TOEddies center at that time, smoothed over 3 days), and the velocity
    orthogonal to the radius -- the azimuthal velocity v_theta, positive anticlockwise -- is
    plotted against distance to the center and depth (bins of R_BIN_KM). The center is
    refined from the SADCP data (see below), because the TOEddies center can be off by
    ~25 km (Subirade et al. 2023), which smears a radial composite. The same refined
    center is used for the LADCP.

Eddy center (CENTER_METHOD)
    "nencioli" (default) -- Nencioli et al. (2008), as applied by Subirade et al. (2023):
        for each time window (one day in radial mode, the whole section in section mode)
        a grid of candidate centers is laid around the SADCP samples and the TOEddies
        center of that window. For every candidate and every depth level of NENC_LAYER the
        velocities are split into radial and tangential components about that candidate,
        and the RMS of the radial velocity is computed (NENC_CRITERION = "rms_radial"; or
        the mean tangential velocity in the eddy's sense of rotation, "mean_tangential").
        The field is averaged over the levels and the center is its minimum. A minimum on
        the edge of the grid is rejected. In radial mode the daily offsets from TOEddies
        are interpolated in time, so the center follows the TOEddies drift between days.
        Diagnostic figure: <eddy>_<section>_center_nencioli.png (cf. Subirade Fig. S2).
    "offset"   -- one constant offset of the TOEddies track over the whole window,
        minimising the radial energy fraction (fit_center), 10 km <= r <= 1.5 R_max.
    "toeddies" -- no refinement.

One SADCP figure and one LADCP figure per section, saved in the eddy's folder in
figures/eddy_tracks/ (found by looking for tracked_eddy_<name>.csv).

Run:  python eddy_velocity_sections.py            (all eddies)
      python eddy_velocity_sections.py CE18       (one eddy)
"""

import os
import sys
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines
import cartopy.crs as ccrs
import cartopy.feature as cfeature

try:
    import cmocean
    CMAP = cmocean.cm.balance
except ImportError:
    CMAP = "RdBu_r"
try:
    import geopandas as gpd
except ImportError:
    gpd = None

# =============================================================================
# Eddies and sections — edit here
# =============================================================================
EDDIES = {
    "CE18": dict(
        name="CE18_202601601", eddy_type=-1,
        sections={
            # straight NE->SW crossing through the center (st. 071 -> 073 -> 074, 13 km from center)
            "crossing_0716": dict(mode="section", t0="2026-07-16T06:00", t1="2026-07-17T16:00"),
            # western flank, 54-89 km from the center (st. 044 -> 069)
            "flank_0710": dict(mode="section", t0="2026-07-10T10:00", t1="2026-07-12T17:00"),
        }),
    "AE11": dict(
        name="AE11_20260514", eddy_type=1,
        sections={
            # stations 005-048: 27 to 186 km from the center, all around the eddy
            "radial_0701_0710": dict(mode="radial", t0="2026-07-01T00:00", t1="2026-07-11T00:00",
                                     r_limit_km=220),
        }),
}

SADCP_AVG = "10min"          # SADCP time averaging (~3 km at 10 kn)
USE_DETIDED = True           # use UVEL_ADCP_CORTIDE / VVEL_ADCP_CORTIDE when present
Z_SPLIT = 150.0              # m: OS150 above, OS38 below
ZMAX_SADCP = 1000.0
ZMAX_LADCP = 1500.0
DS_KM = 3.0                  # bin size along a section (SADCP)
R_BIN_KM = 10.0              # radius bin (radial mode, SADCP)
LADCP_MAX_GAP_KM = 60.0      # do not interpolate between LADCP stations further apart
CENTER_METHOD = "nencioli"   # "nencioli" | "offset" | "toeddies" (see module docstring)
# --- Nencioli et al. (2008) / Subirade et al. (2023)
NENC_CRITERION = "rms_radial"  # "rms_radial": min RMS of v_r (Nencioli)
                               # "mean_tangential": max mean v_theta in the eddy's sense of
                               #   rotation (= Subirade's "minimal tangential" for a NH anticyclone)
NENC_LAYER = (0.0, 90.0)     # m: levels used, cost averaged over them (Subirade: surface-90 m)
NENC_WINDOW = "1D"           # radial mode: one center per window (Subirade: per day)
NENC_RES_KM = 2.0            # candidate grid spacing (Subirade: 1/8 deg ~ 12 km)
NENC_PAD_KM = 60.0           # grid extends this far beyond the samples and the TOEddies center
NENC_NMIN = 12               # min SADCP profiles (10-min averages) in a window / at a level
# --- constant-offset fit (CENTER_METHOD = "offset")
FIT_LAYER = (30.0, 200.0)    # m: layer used to fit the center
FIT_SEARCH_KM = 60.0         # search +- this offset around the TOEddies center
CENTER_SMOOTH_DAYS = 3       # running mean on the daily TOEddies centers (they move on a 1/8° grid)
PROFILE_DEPTHS = [50, 100, 200, 400, 700]   # radial mode: v_theta(r) at these depths

# =============================================================================
# Paths
# =============================================================================
BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
DIR_SADCP = BASE / "data/SADCP"
SADCP_FILES = {"OS150": DIR_SADCP / "WHIRLS_NB150_TOT_osite.nc",   # full cruise (README)
               "OS38": DIR_SADCP / "WHIRLS_BB38_TOT_osite.nc"}
DIR_LADCP = BASE / "data/LADCP/traitement_avec_SADCP"
DIR_SSH = BASE / "data/SSH"
DIR_EDDIES = {"TOEDDIES_better_ADT": BASE / "data/TOEDDIES/TOEDDIES_better_ADT",
              "TOEDDIES": BASE / "data/TOEDDIES/TOEDDIES"}
FNAME_BATHY = BASE / "REP/data/bathymetry/etopo2.nc"
DIR_TRACKS = BASE / "figures/eddy_tracks"
DIR_CODE = BASE / "analysis"           # where ladcp_tools.py lives

R_KM = 6371.0
crs = ccrs.PlateCarree()


# =============================================================================
# Pure geometry (no I/O)
# =============================================================================
def to_xy(lon, lat, lon0, lat0):
    """Local tangent plane, km."""
    x = np.deg2rad(np.asarray(lon, float) - lon0) * R_KM * np.cos(np.deg2rad(lat0))
    y = np.deg2rad(np.asarray(lat, float) - lat0) * R_KM
    return x, y


def fit_line(x, y):
    """Principal axis of the points. Returns origin, along unit vector (oriented in the
    direction of travel) and normal unit vector (to the LEFT of travel)."""
    ok = np.isfinite(x) & np.isfinite(y)
    xm, ym = x[ok].mean(), y[ok].mean()
    w, vec = np.linalg.eigh(np.cov(np.vstack([x[ok] - xm, y[ok] - ym])))
    es = vec[:, np.argmax(w)]
    xi, yi = x[ok][[0, -1]], y[ok][[0, -1]]
    if (xi[1] - xi[0]) * es[0] + (yi[1] - yi[0]) * es[1] < 0:
        es = -es
    en = np.array([-es[1], es[0]])
    return np.array([xm, ym]), es, en


def polar(x, y, u, v):
    """Radius (km), azimuthal (+ anticlockwise) and radial (+ outward) velocity.
    x, y: (n,) ; u, v: (n,) or (n, nz)."""
    r = np.hypot(x, y)
    rr = np.where(r > 0, r, np.nan)
    if np.ndim(u) == 2:
        x, y, rr = x[:, None], y[:, None], rr[:, None]
    vt = (-u * y + v * x) / rr
    vr = (u * x + v * y) / rr
    return r, vt, vr


def fit_center(x, y, u, v, r_max, search=FIT_SEARCH_KM, step=2.0, r_min=10.0, nmin=20):
    """Offset (dx, dy) of the center minimising the fraction of radial kinetic energy."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(u) & np.isfinite(v)
    x, y, u, v = x[ok], y[ok], u[ok], v[ok]
    best = (np.inf, 0.0, 0.0)
    for dx in np.arange(-search, search + step / 2, step):
        for dy in np.arange(-search, search + step / 2, step):
            r, vt, vr = polar(x - dx, y - dy, u, v)
            m = (r >= r_min) & (r <= r_max)
            if m.sum() < nmin:
                continue
            cost = np.sum(vr[m] ** 2) / np.sum(vr[m] ** 2 + vt[m] ** 2)
            if cost < best[0]:
                best = (cost, dx, dy)
    return best[1], best[2], best[0]


def nencioli_cost(x, y, U, V, gx, gy, criterion="rms_radial", rot_sign=1, nmin=NENC_NMIN):
    """Nencioli et al. (2008) cost of every candidate center, at every level.

    x, y   : (n,) sample positions, km ; U, V : (n, nz) velocities, m/s
    gx, gy : (nx,), (ny,) candidate-center coordinates, km
    For a candidate center, each velocity is split into a radial (outward) and a tangential
    (anticlockwise) component, with the angle measured from the candidate to the sample:
        v_r = u cos(a) + v sin(a),   v_t = -u sin(a) + v cos(a)
    Returns cost (nz, ny, nx), to be MINIMISED:
        rms_radial      : sqrt(mean(v_r^2))           (Nencioli: center = min RMS radial)
        mean_tangential : -rot_sign * mean(v_t)       (center = max rotation, rot_sign = +1
                                                       anticlockwise, -1 clockwise)
    Levels with fewer than nmin valid samples are NaN."""
    GX, GY = np.meshgrid(gx, gy)
    DX = np.asarray(x, float)[None, None, :] - GX[..., None]
    DY = np.asarray(y, float)[None, None, :] - GY[..., None]
    R = np.hypot(DX, DY)
    R[R == 0] = np.nan
    ca, sa = DX / R, DY / R
    cost = np.full((U.shape[1],) + GX.shape, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for k in range(U.shape[1]):
            ok = np.isfinite(U[:, k]) & np.isfinite(V[:, k]) & np.isfinite(x) & np.isfinite(y)
            if ok.sum() < nmin:
                continue
            u, v, c, s = U[ok, k], V[ok, k], ca[..., ok], sa[..., ok]
            if criterion == "rms_radial":
                cost[k] = np.sqrt(np.nanmean((u * c + v * s) ** 2, axis=-1))
            elif criterion == "mean_tangential":
                cost[k] = -rot_sign * np.nanmean(-u * s + v * c, axis=-1)
            else:
                raise ValueError(f"unknown NENC_CRITERION {criterion!r}")
    return cost


def nencioli_center(x, y, U, V, criterion="rms_radial", rot_sign=1,
                    res=NENC_RES_KM, pad=NENC_PAD_KM, nmin=NENC_NMIN):
    """Eddy center from one quasi-synoptic set of velocity profiles (Nencioli et al. 2008).

    x, y are km relative to a reference point (here the TOEddies center, at (0, 0)). The
    candidate grid covers the samples and the reference point, padded by `pad`. The cost is
    computed at each level, averaged vertically (Subirade et al. 2023) and minimised.
    Returns a dict: dx, dy (km, center relative to the reference), cost (ny, nx), gx, gy,
    edge (minimum on the grid boundary -> unreliable), spread (km, scatter of the centers
    found level by level), nlev (levels used). None if no level has enough data."""
    ok = np.isfinite(x) & np.isfinite(y)
    gx = np.arange(min(np.nanmin(x[ok]), 0) - pad, max(np.nanmax(x[ok]), 0) + pad + res / 2, res)
    gy = np.arange(min(np.nanmin(y[ok]), 0) - pad, max(np.nanmax(y[ok]), 0) + pad + res / 2, res)
    cost = nencioli_cost(x, y, U, V, gx, gy, criterion, rot_sign, nmin)
    lev = np.isfinite(cost).any(axis=(1, 2))
    if not lev.any():
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        field = np.nanmean(cost[lev], axis=0)
    iy, ix = np.unravel_index(np.nanargmin(field), field.shape)
    edge = ix in (0, len(gx) - 1) or iy in (0, len(gy) - 1)
    # center level by level, as in the original routine (one call per depth)
    pk = np.array([np.unravel_index(np.nanargmin(c), c.shape) for c in cost[lev]])
    spread = float(np.hypot(gx[pk[:, 1]].std(), gy[pk[:, 0]].std())) if len(pk) > 1 else np.nan
    return dict(dx=float(gx[ix]), dy=float(gy[iy]), cost=field, gx=gx, gy=gy, edge=bool(edge),
                spread=spread, nlev=int(lev.sum()))


def bin_columns(coord, vals, edges, nmin=1):
    """Mean of vals (n, nz) in bins of coord. Returns centers, (nbins, nz) mean, counts."""
    idx = np.digitize(coord, edges) - 1
    nb = len(edges) - 1
    out = np.full((nb, vals.shape[1]), np.nan)
    cnt = np.zeros(nb, int)
    for b in range(nb):
        m = idx == b
        cnt[b] = m.sum()
        if cnt[b] >= nmin:
            with np.errstate(all="ignore"):
                out[b] = np.nanmean(vals[m], axis=0)
    return 0.5 * (edges[:-1] + edges[1:]), out, cnt


def interp_columns(coord, vals, grid, max_gap):
    """Linear interpolation of station profiles (n, nz) onto grid, not across gaps."""
    o = np.argsort(coord)
    c, V = np.asarray(coord)[o], vals[o]
    out = np.full((len(grid), V.shape[1]), np.nan)
    for k in range(V.shape[1]):
        m = np.isfinite(V[:, k])
        if m.sum() >= 2:
            cc, vv = c[m], V[m, k]
            col = np.interp(grid, cc, vv, left=np.nan, right=np.nan)
            j = np.clip(np.searchsorted(cc, grid), 1, len(cc) - 1)
            near = np.min(np.abs(grid[:, None] - cc[None, :]), axis=1) <= np.diff(grid).mean() / 2
            col[((cc[j] - cc[j - 1]) > max_gap) & ~near] = np.nan
            out[:, k] = col
        elif m.sum() == 1:
            out[np.argmin(np.abs(grid - c[m][0])), k] = V[m, k][0]
    return out


# =============================================================================
# Data loading
# =============================================================================
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


_SADCP_CACHE = {}


def load_sadcp_raw(path):
    """Full-resolution SADCP as Dataset(time, z) with u, v, lat, lon (QC flag == 1)."""
    if path in _SADCP_CACHE:
        return _SADCP_CACHE[path]
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
    _SADCP_CACHE[path] = out
    return out


def sadcp_window(t0, t1):
    """OS150 above Z_SPLIT + OS38 below, averaged over SADCP_AVG, inside [t0, t1]."""
    parts = []
    for inst, f in SADCP_FILES.items():
        if not f.exists():
            print(f"  ! missing {f}")
            continue
        d = load_sadcp_raw(f).sel(time=slice(t0, t1))
        d = d.resample(time=SADCP_AVG).mean()
        zz = d.z.values
        keep = zz <= Z_SPLIT if inst == "OS150" else (zz > Z_SPLIT) & (zz <= ZMAX_SADCP)
        d = d.isel(z=np.where(keep)[0])
        parts.append(d)
    if not parts:
        raise FileNotFoundError("no SADCP file found")
    base = parts[0]
    out = base
    for p in parts[1:]:
        p = p.reindex(time=base.time)
        out = xr.Dataset({"u": xr.concat([out.u, p.u], "z"), "v": xr.concat([out.v, p.v], "z"),
                          "lat": out.lat.fillna(p.lat), "lon": out.lon.fillna(p.lon)})
    out = out.sortby("z")
    ok = np.isfinite(out.lat.values) & np.isfinite(out.lon.values) & np.isfinite(out.u.values).any(1)
    return out.isel(time=ok)


def load_ladcp(t0, t1):
    """Stations cast in [t0, t1]: list of dicts with station, time, lon, lat, z, u, v."""
    sys.path.insert(0, str(DIR_CODE))
    from ladcp_tools import load_all_stations, read_lad
    st = load_all_stations(str(DIR_LADCP) + "/")
    tcol = "datetime" if "datetime" in st.columns else st.select_dtypes("datetime").columns[0]
    st = st[(st[tcol] >= pd.Timestamp(t0)) & (st[tcol] <= pd.Timestamp(t1))]
    out = []
    for _, s in st.iterrows():
        f = s["file"] if "file" in s and str(s["file"]).endswith(".lad") else \
            os.path.join(str(DIR_LADCP), str(s.station), f"{s.station}.lad")
        try:
            _, df = read_lad(f)
        except Exception as e:
            print(f"  ! LADCP {s.station}: {e}")
            continue
        cols = {c.lower(): c for c in df.columns}
        zc = next((cols[c] for c in ("z", "depth", "dep", "pressure", "p") if c in cols), None)
        uc = next((cols[c] for c in ("u", "u_mean", "uvel") if c in cols), None)
        vc = next((cols[c] for c in ("v", "v_mean", "vvel") if c in cols), None)
        if None in (zc, uc, vc):
            raise KeyError(f"LADCP columns not recognised: {list(df.columns)} "
                           "-> edit the names in load_ladcp()")
        df = df[[zc, uc, vc]].astype(float).dropna()
        out.append(dict(station=str(s.station), time=pd.Timestamp(s[tcol]), lon=float(s.lon),
                        lat=float(s.lat), z=np.abs(df[zc].values), u=df[uc].values,
                        v=df[vc].values))
    print(f"  LADCP: {len(out)} station(s) in the window: {', '.join(o['station'] for o in out)}")
    return out


def ladcp_on_grid(stations, zgrid):
    """(n_station, nz) u and v on a common depth grid."""
    U = np.full((len(stations), len(zgrid)), np.nan)
    V = np.full_like(U, np.nan)
    for i, s in enumerate(stations):
        o = np.argsort(s["z"])
        U[i] = np.interp(zgrid, s["z"][o], s["u"][o], left=np.nan, right=np.nan)
        V[i] = np.interp(zgrid, s["z"][o], s["v"][o], left=np.nan, right=np.nan)
    return U, V


def load_track(name):
    f = next(iter(sorted(DIR_TRACKS.rglob(f"tracked_eddy_{name}.csv"))), None)
    if f is None:
        raise FileNotFoundError(f"tracked_eddy_{name}.csv not found under {DIR_TRACKS}")
    tr = pd.read_csv(f, parse_dates=["date"])
    return tr, f.parent


def center_at(tr, times):
    """TOEddies center at given times: daily centers at 12:00, running mean, linear interp."""
    t = tr.date + pd.Timedelta(hours=12)
    lon = tr.lon.rolling(CENTER_SMOOTH_DAYS, center=True, min_periods=1).mean().values
    lat = tr.lat.rolling(CENTER_SMOOTH_DAYS, center=True, min_periods=1).mean().values
    tq = pd.to_datetime(times).values.astype("datetime64[s]").astype(float)
    tt = t.values.astype("datetime64[s]").astype(float)
    return np.interp(tq, tt, lon), np.interp(tq, tt, lat)


def rmax_at(tr, t0, t1):
    m = (tr.date >= pd.Timestamp(t0).normalize()) & (tr.date <= pd.Timestamp(t1))
    return float(tr.loc[m, "r_maxspeed_km"].mean())


def eddy_contours(tr, day, eddy_type):
    """Max-speed and outermost contours of the tracked eddy on that day (lon/lat arrays)."""
    if gpd is None:
        return None, None
    row = tr.iloc[(tr.date - pd.Timestamp(day).normalize()).abs().argmin()]
    folder = DIR_EDDIES.get(row.source, DIR_EDDIES["TOEDDIES"])
    kind = "anticyclones" if eddy_type == 1 else "cyclones"
    d = row.date
    res = []
    for part in ("maxspeed", "outcontour"):
        xy = None
        for tag in (d.strftime("%Y%m%d"), d.strftime("%Y-%m-%d")):
            for f in sorted(glob.glob(str(folder / f"*{tag}*{part}*"))):
                if (kind == "anticyclones") != ("anticyclon" in Path(f).name.lower()):
                    continue
                g = gpd.read_file(f).geometry
                i = int(row.index_in_file)
                if i < len(g) and g.iloc[i] is not None and not g.iloc[i].is_empty:
                    geom = g.iloc[i]
                    if geom.geom_type.startswith("Multi"):
                        geom = max(geom.geoms, key=lambda q: q.length)
                    c = geom.exterior.coords if geom.geom_type == "Polygon" else geom.coords
                    xy = np.asarray(c)[:, :2]
                break
            if xy is not None:
                break
        res.append(xy)
    return res


def adt_day(day):
    d = pd.Timestamp(day)
    tags = (d.strftime("%Y-%m-%d"), d.strftime("%Y%m%d"))
    f = next((f for f in sorted(glob.glob(str(DIR_SSH / "*.nc")))
              if any(t in Path(f).name for t in tags)), None)
    if f is None:
        return None
    ds = xr.open_dataset(f)
    ds = ds.rename({k: v for k, v in {"latitude": "lat", "longitude": "lon"}.items()
                    if k in ds.variables})
    a = ds["adt"]
    if "time" in a.dims:
        a = a.isel(time=0)
    return a.load()


# =============================================================================
# Plot helpers
# =============================================================================
def vlim_of(*arrs):
    v = np.concatenate([np.abs(np.asarray(a)).ravel() for a in arrs])
    v = v[np.isfinite(v)]
    return max(0.1, np.round(np.nanpercentile(v, 98) * 10) / 10) if v.size else 0.5


def draw_section(ax, xg, z, V, vlim, xlabel, title):
    levels = np.linspace(-vlim, vlim, 21)
    Vm = np.ma.masked_invalid(V.T)
    cf = ax.contourf(xg, z, Vm, levels=levels, cmap=CMAP, extend="both")
    try:
        cs = ax.contour(xg, z, Vm, levels=np.arange(-2, 2.01, 0.1), colors="k", linewidths=0.4)
        ax.clabel(cs, fmt="%.1f", fontsize=6)
    except Exception:
        pass
    ax.set_ylim(z.max() if np.isfinite(z.max()) else 1000, 0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Depth (m)")
    ax.set_title(title, loc="left", fontsize=10)
    return cf


def geo_map(ax, extent, tr, day, eddy_type):
    ax.set_extent(extent, crs=crs)
    ax.add_feature(cfeature.LAND, facecolor="0.6", edgecolor="k", zorder=3)
    gl = ax.gridlines(draw_labels=True, lw=0.3, color="0.5", alpha=0.5)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 7}
    a = adt_day(day)
    if a is not None:
        a = a.sel(lon=slice(extent[0] - 1, extent[1] + 1), lat=slice(extent[2] - 1, extent[3] + 1))
        ax.contour(a.lon, a.lat, a, levels=np.arange(-2, 3, 0.05), colors="0.5",
                   linewidths=0.4, transform=crs, zorder=1)
    ms, oc = eddy_contours(tr, day, eddy_type)
    col = "tab:red" if eddy_type == 1 else "tab:blue"
    if ms is not None:
        ax.plot(ms[:, 0], ms[:, 1], color=col, lw=1.5, transform=crs, zorder=4)
    if oc is not None:
        ax.plot(oc[:, 0], oc[:, 1], color="k", lw=0.8, ls=":", transform=crs, zorder=4)


# =============================================================================
# Eddy center refinement
# =============================================================================
_FITTED = {}


def refine_center(times, lon, lat, z, U, V, tr, cfg, per_window, label, out_dir=None):
    """Center of the eddy relative to the TOEddies track, from SADCP.

    times (n,), lon, lat (n,), z (nz,), U, V (n, nz). Returns dict with t (k,) datetime64 and
    dx, dy (k,) km east/north of the TOEddies center at t; use offset_at() to get the offset
    at any time. per_window: one center per NENC_WINDOW (radial mode) or one for all data
    (section mode). Only CENTER_METHOD == "nencioli" uses per_window."""
    times = pd.DatetimeIndex(pd.to_datetime(times))
    tmid = np.array([times[len(times) // 2]], "datetime64[ns]")
    if CENTER_METHOD == "toeddies":
        return dict(t=tmid, dx=np.zeros(1), dy=np.zeros(1), method="TOEddies center", diags=[])

    clon, clat = center_at(tr, times.values)
    x = np.deg2rad(lon - clon) * R_KM * np.cos(np.deg2rad(clat))
    y = np.deg2rad(lat - clat) * R_KM

    if CENTER_METHOD == "offset":
        lay = (z >= FIT_LAYER[0]) & (z <= FIT_LAYER[1])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            ul, vl = np.nanmean(U[:, lay], 1), np.nanmean(V[:, lay], 1)
        rmax = rmax_at(tr, times[0], times[-1])
        dx, dy, cost = fit_center(x, y, ul, vl, r_max=1.5 * rmax)
        _, vt0, vr0 = polar(x, y, ul, vl)
        m0 = np.hypot(x, y) <= 1.5 * rmax
        c0 = np.nansum(vr0[m0] ** 2) / np.nansum(vr0[m0] ** 2 + vt0[m0] ** 2)
        print(f"  center offset ({dx:+.0f}, {dy:+.0f}) km east/north of TOEddies; "
              f"radial energy fraction {c0:.2f} -> {cost:.2f}")
        return dict(t=tmid, dx=np.array([dx]), dy=np.array([dy]),
                    method="constant offset of TOEddies (min radial energy)", diags=[])

    if CENTER_METHOD != "nencioli":
        raise ValueError(f"unknown CENTER_METHOD {CENTER_METHOD!r}")
    lay = (z >= NENC_LAYER[0]) & (z <= NENC_LAYER[1])
    if not lay.any():
        raise ValueError(f"no SADCP level in NENC_LAYER {NENC_LAYER}")
    Ul, Vl = U[:, lay], V[:, lay]
    # sense of rotation (only used by "mean_tangential"): + anticlockwise
    rot_sign = cfg["eddy_type"] * (1 if np.nanmean(lat) < 0 else -1)
    groups = times.floor(NENC_WINDOW) if per_window else pd.DatetimeIndex([times[0]] * len(times))
    print(f"  Nencioli center ({NENC_CRITERION}, {NENC_LAYER[0]:.0f}-{NENC_LAYER[1]:.0f} m, "
          f"{lay.sum()} levels, grid {NENC_RES_KM:g} km):")
    diags = []
    for g in groups.unique():
        m = np.asarray(groups == g)
        if m.sum() < NENC_NMIN:
            print(f"    {g:%Y-%m-%d}: {m.sum()} profiles < NENC_NMIN, skipped")
            continue
        tk = times[m].mean()
        c_lon, c_lat = (float(a[0]) for a in center_at(tr, [tk]))
        # positions relative to the TOEddies center of THIS window (fixed frame in the window)
        xk = np.deg2rad(lon[m] - c_lon) * R_KM * np.cos(np.deg2rad(c_lat))
        yk = np.deg2rad(lat[m] - c_lat) * R_KM
        res = nencioli_center(xk, yk, Ul[m], Vl[m], NENC_CRITERION, rot_sign)
        if res is None:
            print(f"    {g:%Y-%m-%d}: no level with enough data, skipped")
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            res.update(t=tk, clon=c_lon, clat=c_lat, n=int(m.sum()), x=xk, y=yk,
                       ul=np.nanmean(Ul[m], 1), vl=np.nanmean(Vl[m], 1),
                       lon=c_lon + np.rad2deg(res["dx"] / (R_KM * np.cos(np.deg2rad(c_lat)))),
                       lat=c_lat + np.rad2deg(res["dy"] / R_KM))
        diags.append(res)
        print(f"    {tk:%Y-%m-%d %H:%M}  n={m.sum():4d}  center {res['lon']:.3f}E {res['lat']:.3f}N, "
              f"({res['dx']:+5.0f}, {res['dy']:+5.0f}) km from TOEddies, per-level spread "
              f"{res['spread']:.0f} km{'  ON GRID EDGE -> rejected' if res['edge'] else ''}")
    if out_dir is not None:
        plot_nencioli(diags, cfg, label, out_dir)
    good = [d for d in diags if not d["edge"]]
    if not good:
        print("  ! no valid Nencioli center, falling back to the TOEddies center")
        return dict(t=tmid, dx=np.zeros(1), dy=np.zeros(1), method="TOEddies center", diags=diags)
    return dict(t=np.array([d["t"] for d in good], "datetime64[ns]"),
                dx=np.array([d["dx"] for d in good]), dy=np.array([d["dy"] for d in good]),
                method=f"Nencioli et al. (2008), {NENC_CRITERION}", diags=diags)


def offset_at(cen, times):
    """Center offset (km east, km north of TOEddies) at given times: linear interpolation
    between the window centers, held constant before the first / after the last."""
    tq = pd.to_datetime(times).values.astype("datetime64[s]").astype(float)
    tk = cen["t"].astype("datetime64[s]").astype(float)
    if len(tk) == 1:
        return np.full(tq.shape, cen["dx"][0]), np.full(tq.shape, cen["dy"][0])
    return np.interp(tq, tk, cen["dx"]), np.interp(tq, tk, cen["dy"])


def plot_nencioli(diags, cfg, label, out_dir):
    """One panel per window: cost field, layer-mean velocities, TOEddies and Nencioli centers."""
    if not diags:
        return
    n = len(diags)
    nc = min(4, n)
    nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.3 * nc, 4.1 * nr), squeeze=False)
    clab = ("RMS v$_r$ (m s$^{-1}$)" if NENC_CRITERION == "rms_radial"
            else "$-$mean v$_\\theta$ in the sense of rotation (m s$^{-1}$)")
    for ax, d in zip(axs.flat, diags):
        pc = ax.pcolormesh(d["gx"], d["gy"], d["cost"], cmap="viridis", shading="auto")
        try:
            ax.contour(d["gx"], d["gy"], d["cost"], 12, colors="w", linewidths=0.4)
        except Exception:
            pass
        st = max(1, len(d["x"]) // 25)
        ax.quiver(d["x"][::st], d["y"][::st], d["ul"][::st], d["vl"][::st], angles="xy",
                  scale_units="xy", scale=1 / 40.0, width=0.005, color="w",
                  edgecolor="k", linewidth=0.4)                              # 1 m/s = 40 km
        ax.plot(d["x"], d["y"], "-", color="0.9", lw=0.6)
        ax.plot(0, 0, "x", color="r", ms=8, mew=2)
        ax.plot(d["dx"], d["dy"], "+", color="0.5" if d["edge"] else "m", ms=14, mew=2.5)
        ax.set_aspect("equal")
        ax.set_title(f"{pd.Timestamp(d['t']):%Y-%m-%d %H:%M}  n={d['n']}\n"
                     f"({d['dx']:+.0f}, {d['dy']:+.0f}) km, spread {d['spread']:.0f} km"
                     f"{'  EDGE' if d['edge'] else ''}", fontsize=8)
        ax.tick_params(labelsize=7)
        fig.colorbar(pc, ax=ax, shrink=0.75).set_label(clab, fontsize=7)
    for ax in axs.flat[n:]:
        ax.axis("off")
    for ax in axs[-1]:
        ax.set_xlabel("km east of TOEddies center", fontsize=8)
    for ax in axs[:, 0]:
        ax.set_ylabel("km north of TOEddies center", fontsize=8)
    fig.suptitle(f"{cfg['name']} — {label}: Nencioli et al. (2008) center, {NENC_CRITERION}, "
                 f"{NENC_LAYER[0]:.0f}-{NENC_LAYER[1]:.0f} m mean  "
                 f"(x TOEddies, + Nencioli, arrows: layer-mean SADCP, 1 m/s = 40 km)", fontsize=9)
    fig.tight_layout()
    out = out_dir / f"{cfg['name']}_{label}_center_nencioli.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


def refined_center_for(cfg, key, sad, tr, out_dir, per_window):
    """Cached refine_center on a SADCP window (computed once, reused for the LADCP figure)."""
    ck = (cfg["name"], key)
    if ck not in _FITTED:
        _FITTED[ck] = refine_center(sad.time.values, sad.lon.values, sad.lat.values, sad.z.values,
                                    sad.u.values, sad.v.values, tr, cfg, per_window, key, out_dir)
    return _FITTED[ck]


# =============================================================================
# Section mode
# =============================================================================
def run_section(key, sec, cfg, tr, out_dir, instrument):
    t0, t1 = sec["t0"], sec["t1"]
    tmid = pd.Timestamp(t0) + (pd.Timestamp(t1) - pd.Timestamp(t0)) / 2
    clon, clat = center_at(tr, [tmid])
    clon, clat = float(clon[0]), float(clat[0])
    rmax = rmax_at(tr, t0, t1)

    # ship track defines the line for both instruments (same line -> comparable figures)
    sad = sadcp_window(t0, t1)
    lon0, lat0 = float(np.nanmean(sad.lon)), float(np.nanmean(sad.lat))
    xs, ys = to_xy(sad.lon.values, sad.lat.values, lon0, lat0)
    origin, es, en = fit_line(xs, ys)
    proj_s = lambda x, y: (x - origin[0]) * es[0] + (y - origin[1]) * es[1]
    proj_n = lambda x, y: (x - origin[0]) * en[0] + (y - origin[1]) * en[1]
    cx, cy = to_xy(clon, clat, lon0, lat0)
    s_c, d_c = float(proj_s(cx, cy)), float(proj_n(cx, cy))
    heading = (np.degrees(np.arctan2(es[0], es[1])) + 360) % 360

    # Nencioli center from the SADCP of the whole section (one quasi-synoptic window)
    nen = None
    if CENTER_METHOD == "nencioli":
        cen = refined_center_for(cfg, key, sad, tr, out_dir, per_window=False)
        good = [d for d in cen["diags"] if not d["edge"]]
        if good:
            ncx, ncy = to_xy(good[0]["lon"], good[0]["lat"], lon0, lat0)
            nen = (good[0]["lon"], good[0]["lat"], float(proj_s(ncx, ncy)), float(proj_n(ncx, ncy)))
    # R_max is drawn around the best available center
    s_b, d_b = (nen[2], nen[3]) if nen else (s_c, d_c)

    if instrument == "SADCP":
        s = proj_s(xs, ys)
        vn = sad.u.values * en[0] + sad.v.values * en[1]
        edges = np.arange(np.floor(s.min()), np.ceil(s.max()) + DS_KM, DS_KM)
        xg, V, _ = bin_columns(s, vn, edges)
        z = sad.z.values
        pts_lon, pts_lat, labels = sad.lon.values, sad.lat.values, None
        dist_line = proj_n(xs, ys)
    else:
        stations = load_ladcp(t0, t1)
        if not stations:
            print(f"  no LADCP station in {key}, skipped")
            return
        z = np.arange(0, ZMAX_LADCP + 1, 16.0)
        U, Vv = ladcp_on_grid(stations, z)
        lo = np.array([q["lon"] for q in stations])
        la = np.array([q["lat"] for q in stations])
        x, y = to_xy(lo, la, lon0, lat0)
        s = proj_s(x, y)
        vn = U * en[0] + Vv * en[1]
        xg = np.arange(np.floor(s.min()) - 2, np.ceil(s.max()) + 2.5, 1.0)
        V = interp_columns(s, vn, xg, LADCP_MAX_GAP_KM)
        pts_lon, pts_lat, labels = lo, la, [q["station"] for q in stations]
        dist_line = proj_n(x, y)
        keepz = np.isfinite(V).any(0)
        if keepz.any():
            z, V = z[: np.where(keepz)[0].max() + 1], V[:, : np.where(keepz)[0].max() + 1]

    vlim = cfg.get("vlim") or vlim_of(V)
    fig = plt.figure(figsize=(15, 6))
    gs = GridSpec(1, 2, width_ratios=[1, 2.3], wspace=0.15, figure=fig)
    axm = fig.add_subplot(gs[0], projection=crs)
    lon_all = np.r_[sad.lon.values, clon]
    lat_all = np.r_[sad.lat.values, clat]
    pad = max(0.5, rmax / 111 * 1.3)
    extent = [np.nanmin(lon_all) - pad, np.nanmax(lon_all) + pad,
              np.nanmin(lat_all) - pad, np.nanmax(lat_all) + pad]
    geo_map(axm, extent, tr, tmid, cfg["eddy_type"])
    axm.plot(sad.lon, sad.lat, "-", color="0.3", lw=1, transform=crs, zorder=5)
    axm.scatter(pts_lon, pts_lat, s=8 if labels is None else 30,
                c="k" if labels is None else "lime", edgecolors="k", lw=0.3,
                transform=crs, zorder=6, marker="." if labels is None else "^")
    if labels:
        for lo_, la_, lab in zip(pts_lon, pts_lat, labels):
            axm.text(lo_ + 0.03, la_ + 0.03, lab, fontsize=6, transform=crs, zorder=7)
    # fitted line and direction of positive orthogonal velocity
    smin, smax = np.nanmin(proj_s(xs, ys)), np.nanmax(proj_s(xs, ys))
    lx = origin[0] + np.array([smin, smax]) * es[0]
    ly = origin[1] + np.array([smin, smax]) * es[1]
    llon = lon0 + np.rad2deg(lx / (R_KM * np.cos(np.deg2rad(lat0))))
    llat = lat0 + np.rad2deg(ly / R_KM)
    tf = crs._as_mpl_transform(axm)
    axm.annotate("", xy=(llon[1], llat[1]), xytext=(llon[0], llat[0]), xycoords=tf,
                 textcoords=tf, arrowprops=dict(arrowstyle="->", color="m", lw=1.5), zorder=8)
    mx, my = llon.mean(), llat.mean()
    nlen = (smax - smin) * 0.15
    axm.annotate("+v$_\\perp$", xy=(mx + np.rad2deg(en[0] * nlen / (R_KM * np.cos(np.deg2rad(lat0)))),
                                     my + np.rad2deg(en[1] * nlen / R_KM)),
                 xytext=(mx, my), xycoords=tf, textcoords=tf,
                 arrowprops=dict(arrowstyle="->", color="m"), color="m", fontsize=8, zorder=8)
    axm.plot(clon, clat, "x", color="k", ms=8, mew=2, transform=crs, zorder=9)
    handles = [
        mlines.Line2D([], [], color="tab:red" if cfg["eddy_type"] == 1 else "tab:blue",
                      label=f"Max-speed contour {tmid:%m-%d}"),
        mlines.Line2D([], [], color="k", ls=":", label="Outermost contour"),
        mlines.Line2D([], [], color="0.3", label="Ship track"),
        mlines.Line2D([], [], color="k", marker="x", ls="", label="TOEddies center"),
        mlines.Line2D([], [], color="m", label="Section axis")]
    if nen:
        axm.plot(nen[0], nen[1], "+", color="g", ms=12, mew=2.5, transform=crs, zorder=9)
        handles.append(mlines.Line2D([], [], color="g", marker="+", ls="", mew=2,
                                     label="Nencioli center (SADCP)"))
    axm.set_title(f"{cfg['name']} — {key}\nline heading {heading:.0f}°, "
                  f"{'Nencioli' if nen else 'TOEddies'} center {abs(d_b):.0f} km from line",
                  loc="left", fontsize=9)
    axm.legend(handles=handles, fontsize=6, loc="lower left")

    ax = fig.add_subplot(gs[1])
    cf = draw_section(ax, xg, z, V, vlim, "Distance along section (km)",
                      f"{instrument}: velocity orthogonal to the section "
                      f"(+ to the left of travel, i.e. towards {(heading - 90) % 360:.0f}°)")
    ax.axvline(s_c, color="k", ls="--", lw=1)
    ax.text(s_c, 0, " TOEddies", va="bottom", fontsize=7)
    if nen:
        ax.axvline(nen[2], color="g", ls="-.", lw=1.2)
        ax.text(nen[2], 0, " Nencioli", va="bottom", fontsize=7, color="g")
    if np.isfinite(rmax) and rmax > abs(d_b):
        half = np.sqrt(rmax ** 2 - d_b ** 2)
        for xx in (s_b - half, s_b + half):
            ax.axvline(xx, color="tab:red" if cfg["eddy_type"] == 1 else "tab:blue", lw=1.5)
        ax.text(s_b + half, 0, " R$_{max}$", va="bottom", fontsize=7)
    if labels:
        for sv, lab in zip(s, labels):
            ax.annotate(lab, (sv, 0), xytext=(0, 2), textcoords="offset points", fontsize=6,
                        ha="center", va="bottom", rotation=90, annotation_clip=False)
            ax.axvline(sv, color="0.3", lw=0.4, ls=":")
    cb = fig.colorbar(cf, ax=ax, pad=0.01)
    cb.set_label("v$_\\perp$ (m s$^{-1}$)")
    fig.text(0.99, 0.01, f"{pd.Timestamp(t0):%Y-%m-%d %H:%M} → {pd.Timestamp(t1):%Y-%m-%d %H:%M}"
             f" | R$_{{max}}$ (TOEddies) {rmax:.0f} km | samples within "
             f"{np.nanmax(np.abs(dist_line)):.0f} km of the section axis", ha="right", fontsize=7)
    out = out_dir / f"{cfg['name']}_{key}_{instrument}_orthogonal_velocity_nencioli.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


# =============================================================================
# Radial mode
# =============================================================================
_SADCP_WIN = {}


def run_radial(key, sec, cfg, tr, out_dir, instrument):
    t0, t1 = sec["t0"], sec["t1"]
    rlim = sec.get("r_limit_km", 200.0)
    rmax = rmax_at(tr, t0, t1)
    tmid = pd.Timestamp(t0) + (pd.Timestamp(t1) - pd.Timestamp(t0)) / 2

    # center refinement from SADCP (done once, reused for the LADCP figure)
    if (cfg["name"], key) not in _SADCP_WIN:
        _SADCP_WIN[(cfg["name"], key)] = sadcp_window(t0, t1)
    sad = _SADCP_WIN[(cfg["name"], key)]
    cen = refined_center_for(cfg, key, sad, tr, out_dir, per_window=True)

    def eddy_frame(times, lon, lat):
        """Position relative to the refined center: TOEddies center at that time (smoothed)
        + refined offset at that time (interpolated between the window centers)."""
        clon, clat = center_at(tr, times)
        x = np.deg2rad(lon - clon) * R_KM * np.cos(np.deg2rad(clat))
        y = np.deg2rad(lat - clat) * R_KM
        ox, oy = offset_at(cen, times)
        return x - ox, y - oy

    if instrument == "SADCP":
        x, y = eddy_frame(sad.time.values, sad.lon.values, sad.lat.values)
        r, vt, _ = polar(x, y, sad.u.values, sad.v.values)
        edges = np.arange(0, rlim + R_BIN_KM, R_BIN_KM)
        rg, V, cnt = bin_columns(r, vt, edges, nmin=3)
        z = sad.z.values
        labels, ptimes = None, sad.time.values
    else:
        stations = load_ladcp(t0, t1)
        if not stations:
            print(f"  no LADCP station in {key}, skipped")
            return
        z = np.arange(0, ZMAX_LADCP + 1, 16.0)
        U, Vv = ladcp_on_grid(stations, z)
        times = np.array([q["time"] for q in stations], "datetime64[ns]")
        lo = np.array([q["lon"] for q in stations])
        la = np.array([q["lat"] for q in stations])
        x, y = eddy_frame(times, lo, la)
        r, vt, _ = polar(x, y, U, Vv)
        rg = np.arange(0, rlim + 1, 2.0)
        # average stations at nearly the same radius (repeat stations), then interpolate
        rb, Vb, cnt = bin_columns(r, vt, np.arange(0, rlim + 5.1, 5.0))
        V = interp_columns(rb[cnt > 0], Vb[cnt > 0], rg, LADCP_MAX_GAP_KM)
        labels, ptimes = [q["station"] for q in stations], times
        keepz = np.isfinite(V).any(0)
        if keepz.any():
            z, V = z[: np.where(keepz)[0].max() + 1], V[:, : np.where(keepz)[0].max() + 1]

    vlim = cfg.get("vlim") or vlim_of(V)
    sign = "anticlockwise" if cfg["eddy_type"] == 1 else "clockwise"
    fig = plt.figure(figsize=(16, 9))
    gs = GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[1, 1.4], hspace=0.3, wspace=0.18,
                  figure=fig)

    # (a) coverage in the eddy frame
    axf = fig.add_subplot(gs[:, 0])
    tnum = (pd.to_datetime(ptimes) - pd.Timestamp(t0)) / pd.Timedelta(days=1)
    sc = axf.scatter(x, y, c=tnum, cmap="viridis", s=6 if labels is None else 40,
                     marker="." if labels is None else "^", edgecolors="none" if labels is None else "k")
    if labels:
        for xx, yy, lab in zip(x, y, labels):
            axf.text(xx + 3, yy + 3, lab, fontsize=6)
    th = np.linspace(0, 2 * np.pi, 200)
    col = "tab:red" if cfg["eddy_type"] == 1 else "tab:blue"
    axf.plot(rmax * np.cos(th), rmax * np.sin(th), color=col, lw=1.5, label=f"R$_{{max}}$ {rmax:.0f} km")
    for rr in np.arange(50, rlim + 1, 50):
        axf.plot(rr * np.cos(th), rr * np.sin(th), color="0.7", lw=0.5)
    axf.plot(0, 0, "+", color="k", ms=12, mew=2, label="refined center")
    # TOEddies center seen from the refined center, one marker per window
    axf.plot(-cen["dx"], -cen["dy"], "x", color="0.4", ms=8, mew=1.5, ls="",
             label="TOEddies center" + (" (per window)" if len(cen["dx"]) > 1 else ""))
    axf.set_aspect("equal")
    axf.set_xlim(-rlim, rlim)
    axf.set_ylim(-rlim, rlim)
    axf.set_xlabel("East of eddy center (km)")
    axf.set_ylabel("North of eddy center (km)")
    axf.legend(fontsize=7, loc="lower left")
    axf.set_title(f"(a) {instrument} positions in the eddy frame", loc="left", fontsize=10)
    cbt = fig.colorbar(sc, ax=axf, orientation="horizontal", pad=0.08, shrink=0.8)
    cbt.set_label(f"days since {pd.Timestamp(t0):%Y-%m-%d}")

    # (b) v_theta(r) at a few depths
    axp = fig.add_subplot(gs[0, 1])
    cmap = plt.get_cmap("viridis")
    for k, zz in enumerate(PROFILE_DEPTHS):
        iz = np.argmin(np.abs(z - zz))
        if abs(z[iz] - zz) > 30:
            continue
        axp.plot(rg, V[:, iz], "-o", ms=2.5, lw=1.2, color=cmap(k / max(1, len(PROFILE_DEPTHS) - 1)),
                 label=f"{z[iz]:.0f} m")
    axp.axhline(0, color="k", lw=0.6)
    axp.axvline(rmax, color=col, lw=1.2)
    axp.set_xlim(0, rlim)
    axp.set_ylabel("v$_\\theta$ (m s$^{-1}$)")
    axp.legend(fontsize=7, ncol=len(PROFILE_DEPTHS))
    axp.grid(alpha=0.3)
    axp.set_title(f"(b) Azimuthal velocity vs radius (+ anticlockwise; this eddy turns {sign})",
                  loc="left", fontsize=10)

    # (c) v_theta(r, z)
    ax = fig.add_subplot(gs[1, 1], sharex=axp)
    cf = draw_section(ax, rg, z, V, vlim, "Distance to eddy center (km)",
                      f"(c) {instrument}: velocity orthogonal to the radius, v$_\\theta$(r, z)")
    ax.axvline(rmax, color=col, lw=1.5)
    if labels:
        for rv, lab in zip(r, labels):
            ax.annotate(lab, (rv, 0), xytext=(0, 2), textcoords="offset points", fontsize=6,
                        ha="center", va="bottom", rotation=90, annotation_clip=False)
            ax.axvline(rv, color="0.3", lw=0.4, ls=":")
    cb = fig.colorbar(cf, ax=ax, pad=0.01)
    cb.set_label("v$_\\theta$ (m s$^{-1}$)")
    fig.suptitle(f"{cfg['name']} — {key}: {pd.Timestamp(t0):%Y-%m-%d} → {pd.Timestamp(t1):%Y-%m-%d}, "
                 f"eddy frame (TOEddies center smoothed over {CENTER_SMOOTH_DAYS} d, refined: "
                 f"{cen['method']}, {len(cen['dx'])} center(s), mean offset "
                 f"({np.mean(cen['dx']):+.0f}, {np.mean(cen['dy']):+.0f}) km)", fontsize=10)
    out = out_dir / f"{cfg['name']}_{key}_{instrument}_orthogonal_velocity_nencioli.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


# =============================================================================
def main(which):
    for key in which:
        cfg = EDDIES[key]
        tr, out_dir = load_track(cfg["name"])
        print(f"\n=== {key} ({cfg['name']}) -> {out_dir}")
        for skey, sec in cfg["sections"].items():
            print(f"-- {skey}: {sec['mode']} {sec['t0']} -> {sec['t1']}")
            run = run_section if sec["mode"] == "section" else run_radial
            for instrument in ("SADCP", "LADCP"):
                try:
                    run(skey, sec, cfg, tr, out_dir, instrument)
                except Exception:
                    import traceback
                    print(f"  ! {instrument} figure for {skey} failed:")
                    traceback.print_exc()


if __name__ == "__main__":
    main(sys.argv[1:] or list(EDDIES))