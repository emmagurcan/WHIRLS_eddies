"""
WHIRLS_adt_isolines_ctd.py
==========================

Static map of the ADT field for a chosen TARGET DATE (default 2026-07-13)
with its ADT ISOLINES every 0.1 m, the geostrophic current vectors, the
TOEddies contours, and the CTD station positions:

    Marion Dufresne (MD) CTDs   small ORANGE bullets   <- data/CTD/ctd_metadata.nc
    SA Agulhas (AG) CTDs        small VIOLET bullets   <- optional .xlsx (skipped if absent)

Adapted to Emma's WHIRLS folders (paths in the CONFIG block):

    data/SSH/                 daily CMEMS files (adt, ugos, vgos); the file whose
                              date (YYYY-MM-DD or YYYYMMDD in the name) is NEAREST
                              the target date is used
    data/TOEDDIES/...         TOEddies shapefiles / geojson; better-ADT product
                              first, standard TOEddies as fallback
    data/bathymetry/etopo2.nc bathymetry (variable 'topo')
    data/CTD/ctd_metadata.nc  MD CTD positions (latitude / longitude per cast)

Figures (saved in figures/):
    <date>_WHIRLS_adt_isolines_ctd.png          full region
    <date>_WHIRLS_adt_isolines_ctd_zoom.png     zoom

Usage:
    python WHIRLS_adt_isolines_ctd.py                 target date 2026-07-13
    python WHIRLS_adt_isolines_ctd.py 2026-07-15      another target date
"""
from cmcrameri import cm as cmc
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe

import cartopy.crs as ccrs
import cartopy.feature as cfeature


# ===========================================================================
#  CONFIG
# ===========================================================================
TARGET_DATE = "2026-07-13"

BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
DIR_SSH = BASE / "data/SSH"
DIR_EDDIES = {"TOEDDIES_better_ADT": BASE / "data/TOEDDIES/TOEDDIES_better_ADT",
              "TOEDDIES": BASE / "data/TOEDDIES/TOEDDIES"}
EDDY_PRODUCT = "TOEDDIES_better_ADT"        # tried first; the other one is the fallback
FNAME_BATHY = BASE / "data/bathymetry/etopo2.nc"
DIR_CTD = BASE / "data/CTD"
MD_CTD_NC = DIR_CTD / "ctd_metadata.nc"
OUT_DIR = BASE / "figures"

# SA Agulhas CTDs: give the .xlsx path if you have it. None = look for an
# 'Agu*.xlsx' file anywhere under data/CTD; if none is found they are skipped.
AGU_CTD_XLSX = None
AGU_CTD_SHEET = "Profile"

LON_MIN, LON_MAX = 6.0, 22.0
LAT_MIN, LAT_MAX = -42.0, -30.0
ZOOM_EXTENT = (8.0, 16.0, -40.0, -34.0)

# CTD bullet styling: MD orange / Agulhas violet, thin dark edge.
MD_CTD_COLOR = "#FF8C00"
AGU_CTD_COLOR = "#9400D3"
CTD_MS = 3.5
CTD_MEW = 0.3
CTD_MEC = "k"
CTD_ZORDER = 8

BATHY_LEVELS = [100, 200, 500, 1000, 2000, 3000, 4000, 5000]
BATHY_COLOR = "0.45"
BATHY_LW = 0.35
BATHY_LW_DEEP_FACTOR = 1.25
BATHY_LS = "solid"
BATHY_ALPHA = 0.7
BATHY_ZORDER = 2.5
BATHY_LABELS = False

# Eddy contours: bold black, white halo; cyclones solid, anticyclones dashed
EDDY_COLOR = "k"
EDDY_LW = 1.3
EDDY_ZORDER = 5
EDDY_HALO = True
EDDY_HALO_COLOR = "white"
EDDY_HALO_LW = 2.4
EDDY_HALO_ALPHA = 0.85
EDDY_CENTER_MS = 6

ADT_VMIN, ADT_VMAX = -0.7, 1.7
ADT_CMAP_NAME = "vik"       # Crameri 'vik' if cmcrameri is installed, else cmocean balance

# ADT isolines every 0.1 m
ADT_ISO_STEP = 0.1
ADT_ISO_COLOR = "0.15"
ADT_ISO_LW = 0.5
ADT_ISO_ALPHA = 0.85
ADT_ISO_ZORDER = 3.0
ADT_ISO_LABELS = True
ADT_ISO_LABEL_FS = 5

QUIVER_STEP = 2
QUIVER_STEP_ZOOM = 1
QUIVER_SCALE = 110
QUIVER_SCALE_ZOOM = 45
QUIVER_WIDTH = 0.0007
QUIVER_WIDTH_ZOOM = 0.0012
QUIVER_HEADW = 3.0
QUIVER_HEADL = 3.5
QUIVER_ALPHA = 0.9
QUIVER_COLOR = "k"
QUIVER_KEY_SPEED = 1.0

TICK_MAJOR_FULL = 2.0
TICK_MAJOR_ZOOM = 1.0
N_MINOR = 4
GRID_COLOR = "white"
GRID_ALPHA = 0.45
GRID_LW = 0.4
FIG_DPI = 300
FIG_SIZE_FULL = (13, 9)
FIG_SIZE_ZOOM = (12, 9)
# ===========================================================================

# ---------------------------------------------------------------------------
#  Dates in file names
# ---------------------------------------------------------------------------
_DATE_RE = re.compile(r"(20\d{2})-?(\d{2})-?(\d{2})")


def file_date(path):
    """First YYYY-MM-DD or YYYYMMDD in the file name -> datetime (None if none)."""
    m = _DATE_RE.search(Path(path).name)
    if not m:
        return None
    try:
        return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    except ValueError:
        return None


def nearest_date(dates, target):
    dates = [d for d in dates if d is not None]
    if not dates:
        return None
    return min(dates, key=lambda d: abs((d - target).days))


# ---------------------------------------------------------------------------
#  SSH
# ---------------------------------------------------------------------------
def _subset_bounds(ds, name, lo, hi):
    """Slice a coordinate to [lo, hi] whatever its storage order."""
    vals = ds[name].values
    if vals.size > 1 and vals[0] > vals[-1]:
        return {name: slice(hi, lo)}
    return {name: slice(lo, hi)}


def ssh_for_date(target):
    """SSH file in DIR_SSH NEAREST `target` -> (LON, LAT, adt, ugos, vgos, 'YYYY-MM-DD')."""
    try:
        tgt = datetime.strptime(target, "%Y-%m-%d")
    except ValueError:
        raise SystemExit(f"TARGET_DATE '{target}' is not YYYY-MM-DD.")

    files = [f for f in sorted(DIR_SSH.glob("*.nc")) if file_date(f) is not None]
    if not files:
        raise SystemExit(f"No dated .nc files in {DIR_SSH}")
    path = min(files, key=lambda f: abs((file_date(f) - tgt).days))
    d = file_date(path)
    stem = d.strftime("%Y-%m-%d")
    gap = abs((d - tgt).days)
    if gap == 0:
        print(f"SSH file for {target}: {path}")
    else:
        print(f"[warn] no SSH file for {target}; using the nearest, {stem} "
              f"({gap} day{'s' if gap != 1 else ''} away): {path}")

    ds = xr.open_dataset(path)
    sel = {}
    sel.update(_subset_bounds(ds, "longitude", LON_MIN, LON_MAX))
    sel.update(_subset_bounds(ds, "latitude", LAT_MIN, LAT_MAX))
    sub = ds.sel(**sel)

    nlat, nlon = sub.sizes.get("latitude", 0), sub.sizes.get("longitude", 0)
    if nlat < 2 or nlon < 2:
        full_lat, full_lon = ds["latitude"].values, ds["longitude"].values
        ds.close()
        raise SystemExit(
            f"\nSSH subset is empty or degenerate: {nlat} lat x {nlon} lon.\n"
            f"  requested window : lon {LON_MIN}..{LON_MAX}, lat {LAT_MIN}..{LAT_MAX}\n"
            f"  file covers      : lon {full_lon.min():.3f}..{full_lon.max():.3f}, "
            f"lat {full_lat.min():.3f}..{full_lat.max():.3f}\n"
            f"Adjust LON_MIN/LON_MAX/LAT_MIN/LAT_MAX to overlap the file.")

    LON, LAT = sub["longitude"].values, sub["latitude"].values
    adt = np.atleast_2d(np.squeeze(sub["adt"].values))
    u = np.atleast_2d(np.squeeze(sub["ugos"].values))
    v = np.atleast_2d(np.squeeze(sub["vgos"].values))
    ds.close()

    if LAT.size > 1 and LAT[0] > LAT[-1]:
        LAT = LAT[::-1]
        adt, u, v = adt[::-1, :], u[::-1, :], v[::-1, :]
    if LON.size > 1 and LON[0] > LON[-1]:
        LON = LON[::-1]
        adt, u, v = adt[:, ::-1], u[:, ::-1], v[:, ::-1]

    print(f"  grid {LAT.size} lat x {LON.size} lon; lon {LON.min():.3f}..{LON.max():.3f}, "
          f"lat {LAT.min():.3f}..{LAT.max():.3f}")
    return LON, LAT, adt, u, v, stem


# ---------------------------------------------------------------------------
#  TOEddies (shapefiles or geojson, read with geopandas)
# ---------------------------------------------------------------------------
def _geom_coords(geom):
    """Polygon -> exterior ring; MultiPolygon -> first polygon (as the original);
    Point -> (x, y)."""
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type == "Polygon":
        return list(geom.exterior.coords)
    if geom.geom_type == "MultiPolygon":
        return list(geom.geoms[0].exterior.coords)
    if geom.geom_type == "Point":
        return (geom.x, geom.y)
    if geom.geom_type in ("LineString", "LinearRing"):
        return list(geom.coords)
    return None


def eddies_for_date(ssh_date):
    """Eddy geometry for the plotted SSH date (nearest available date, with warning)."""
    import geopandas as gpd

    tgt = datetime.strptime(ssh_date, "%Y-%m-%d")
    order = [EDDY_PRODUCT] + [k for k in DIR_EDDIES if k != EDDY_PRODUCT]
    for product in order:
        d = DIR_EDDIES[product]
        if not d.exists():
            continue
        files = [f for f in d.iterdir()
                 if f.suffix.lower() in (".shp", ".geojson", ".json") and file_date(f)]
        best = nearest_date([file_date(f) for f in files], tgt)
        if best is None:
            continue
        gap = abs((best - tgt).days)
        if gap > 0 and product == EDDY_PRODUCT:
            # the exact day may exist in the fallback product
            other = [DIR_EDDIES[k] for k in order[1:] if DIR_EDDIES[k].exists()]
            if any(file_date(f) == tgt for o in other for f in o.iterdir()):
                print(f"[info] no {product} file for {ssh_date}; trying the fallback product")
                continue
        files = [f for f in files if file_date(f) == best]
        if gap:
            print(f"[warn] no eddy file for the SSH date {ssh_date}; using "
                  f"{best:%Y-%m-%d} ({gap} day{'s' if gap != 1 else ''} away)")
        print(f"Eddies: {product}, {best:%Y-%m-%d} ({len(files)} files)")

        out = {k: [] for k in ("ac_out", "cy_out", "ac_ms", "cy_ms", "ac_cen", "cy_cen")}
        for f in files:
            name = f.name.lower()
            kind = "ac" if "anti" in name else "cy"
            part = ("out" if "outcontour" in name else "ms" if "maxspeed" in name
                    else "cen" if "center" in name else None)
            if part is None:
                continue
            for geom in gpd.read_file(f).geometry:
                c = _geom_coords(geom)
                if c is not None:
                    out[f"{kind}_{part}"].append(c)
        return out

    print("[warn] no TOEddies files found; eddies skipped.")
    return {}


# ---------------------------------------------------------------------------
#  CTD positions
# ---------------------------------------------------------------------------
def _find_position_var(ds, short, std_name, units):
    """Name of the lat (or lon) variable of the cast START position.
    Matches the name, the standard_name or the units, and prefers names without
    'end' / 'bottom' (ctd_metadata.nc also stores end and bottom positions)."""
    cands = []
    for n, v in ds.variables.items():
        nl = n.lower()
        sn = str(v.attrs.get("standard_name", "")).lower()
        un = str(v.attrs.get("units", "")).lower()
        if nl in (short, std_name) or nl.startswith(short) or sn == std_name or un == units:
            cands.append(n)
    if not cands:
        return None

    def rank(n):
        nl = n.lower()
        return (("end" in nl) or ("bottom" in nl),       # start position first
                nl not in (short, std_name),             # exact names first
                "start" not in nl)
    return sorted(cands, key=rank)[0]


def load_ctd_positions_nc(path, label="MD CTD"):
    """Cast positions from ctd_metadata.nc (variables latitude / longitude)."""
    if not Path(path).exists():
        print(f"[warn] {label}: file not found: {path}; no positions drawn.")
        return np.array([]), np.array([])
    ds = xr.open_dataset(path)
    lat = _find_position_var(ds, "lat", "latitude", "degrees_north")
    lon = _find_position_var(ds, "lon", "longitude", "degrees_east")
    if lat is None or lon is None:
        ds.close()
        raise SystemExit(f"{label}: no latitude/longitude variable found in {path}.\n"
                         f"Variables: {list(ds.variables)}")
    print(f"{label}: positions from variables '{lat}' / '{lon}'")
    lats = np.asarray(ds[lat].values, float).ravel()
    lons = np.asarray(ds[lon].values, float).ravel()
    span = ""
    if "time_start" in ds.variables:
        t = ds["time_start"].values
        if np.issubdtype(t.dtype, np.datetime64):
            t = t[~np.isnat(t)]
            if t.size:
                span = (f"; dates {np.datetime_as_string(t.min(), 'D')} .. "
                        f"{np.datetime_as_string(t.max(), 'D')}")
    ds.close()
    ok = np.isfinite(lats) & np.isfinite(lons)
    print(f"{label}: {ok.sum()} CTD positions from {path}{span}")
    return lons[ok], lats[ok]


def _excel_serial_to_date(x):
    try:
        return datetime(1899, 12, 30) + timedelta(days=float(x))
    except (TypeError, ValueError):
        return None


def find_agulhas_xlsx():
    if AGU_CTD_XLSX:
        return Path(AGU_CTD_XLSX)
    hits = sorted(p for p in DIR_CTD.rglob("*.xlsx") if "agu" in p.name.lower())
    if hits:
        print(f"[info] SA Agulhas CTD file found: {hits[0]}")
        return hits[0]
    print("[info] no SA Agulhas CTD file (Agu*.xlsx) under data/CTD; "
          "Agulhas CTDs not drawn. Set AGU_CTD_XLSX if you have it elsewhere.")
    return None


def load_ctd_positions_xlsx(path, sheet=None, dedupe_by_station=False, label="CTD"):
    """
    CTD positions from an .xlsx file -> (lons, lats). Header row found
    automatically (first row with a cell starting with 'lat'); columns matched
    case-insensitively. With dedupe_by_station, consecutive rows of the same
    station (one row per depth bin) collapse to one position.
    """
    if path is None:
        return np.array([]), np.array([])
    if not Path(path).exists():
        print(f"[warn] {label}: file not found: {path}; no positions drawn.")
        return np.array([]), np.array([])
    try:
        from openpyxl import load_workbook
    except ImportError:
        print(f"[warn] {label}: reading .xlsx needs openpyxl "
              "(python -m pip install openpyxl); skipped.")
        return np.array([]), np.array([])

    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet] if (sheet and sheet in wb.sheetnames) else wb[wb.sheetnames[0]]
    if sheet and sheet not in wb.sheetnames:
        print(f"[warn] {label}: sheet '{sheet}' not in {path}; using '{ws.title}'. "
              f"Sheets: {', '.join(wb.sheetnames)}")
    rows = ws.iter_rows(values_only=True)

    header, i_lat = None, None
    for r in rows:
        cells = ["" if c is None else str(c).strip().lower() for c in r]
        for i, c in enumerate(cells):
            if c.startswith("lat"):
                header, i_lat = cells, i
                break
        if header is not None:
            break
    if header is None:
        wb.close()
        print(f"[warn] {label}: no latitude column found in {path}; skipped.")
        return np.array([]), np.array([])

    def find(*prefixes):
        for i, c in enumerate(header):
            if any(c.startswith(p) for p in prefixes):
                return i
        return None

    i_lon = find("lon", "long")
    i_stn = find("station", "cast", "point")
    i_date = find("date", "mm/dd")
    if i_lon is None:
        wb.close()
        print(f"[warn] {label}: no longitude column found in {path}; skipped.")
        return np.array([]), np.array([])

    lats, lons, dates, last_key = [], [], [], object()
    for r in rows:
        if r is None or max(i_lat, i_lon) >= len(r):
            continue
        try:
            la, lo = float(r[i_lat]), float(r[i_lon])
        except (TypeError, ValueError):
            continue
        key = (r[i_stn] if (i_stn is not None and i_stn < len(r))
               else (round(la, 4), round(lo, 4)))
        if dedupe_by_station and key == last_key:
            continue
        last_key = key
        lats.append(la)
        lons.append(lo)
        if i_date is not None and i_date < len(r):
            d = _excel_serial_to_date(r[i_date])
            if d:
                dates.append(d)
    wb.close()

    lats, lons = np.array(lats, float), np.array(lons, float)
    if lats.size and "[s]" in header[i_lat] and (lats > 0).all():
        print(f"[warn] {label}: latitude header says [S] but values are positive; flipping sign.")
        lats = -lats
    if lons.size and "[w]" in header[i_lon] and (lons > 0).all():
        print(f"[warn] {label}: longitude header says [W] but values are positive; flipping sign.")
        lons = -lons
    span = f"; dates {min(dates):%Y-%m-%d} .. {max(dates):%Y-%m-%d}" if dates else ""
    print(f"{label}: {lats.size} CTD positions from {path} (sheet '{ws.title}'){span}")
    return lons, lats


# ---------------------------------------------------------------------------
#  Drawing
# ---------------------------------------------------------------------------
def draw_adt_isolines(ax, LON, LAT, adt):
    if not np.isfinite(adt).any():
        print("[warn] ADT field is all-NaN; isolines skipped.")
        return None
    lo = np.floor(np.nanmin(adt) / ADT_ISO_STEP) * ADT_ISO_STEP
    hi = np.ceil(np.nanmax(adt) / ADT_ISO_STEP) * ADT_ISO_STEP
    levels = np.arange(lo, hi + ADT_ISO_STEP / 2.0, ADT_ISO_STEP)
    cs = ax.contour(LON, LAT, adt, levels=levels, colors=ADT_ISO_COLOR,
                    linewidths=ADT_ISO_LW, linestyles="solid",
                    transform=ccrs.PlateCarree(), zorder=ADT_ISO_ZORDER, alpha=ADT_ISO_ALPHA)
    if ADT_ISO_LABELS:
        ax.clabel(cs, fmt="%.1f", fontsize=ADT_ISO_LABEL_FS, inline=True, inline_spacing=2)
    print(f"ADT isolines: {levels.size} levels, {levels[0]:.1f} .. {levels[-1]:.1f} m, "
          f"step {ADT_ISO_STEP:g} m")
    return cs


def draw_bathymetry(ax, extent):
    if not FNAME_BATHY.exists():
        print(f"[warn] bathymetry not found: {FNAME_BATHY}")
        return
    lo0, lo1, la0, la1 = extent
    ds = xr.open_dataset(FNAME_BATHY)
    dv = next((v for v in ("topo", "elevation", "z", "Band1", "depth") if v in ds.variables), None)
    if dv is None:
        ds.close()
        print(f"[warn] no elevation variable in {FNAME_BATHY}; bathymetry skipped.")
        return
    latn = "lat" if "lat" in ds.coords else "latitude"
    lonn = "lon" if "lon" in ds.coords else "longitude"
    lat, lon = ds[latn].values, ds[lonn].values
    la_ok = (lat >= la0 - 1) & (lat <= la1 + 1)
    lo_ok = (lon >= lo0 - 1) & (lon <= lo1 + 1)
    elev = ds[dv].isel({latn: np.flatnonzero(la_ok), lonn: np.flatnonzero(lo_ok)})
    elev = np.squeeze(elev.transpose(latn, lonn).values).astype(float)
    ds.close()
    lat, lon = lat[la_ok], lon[lo_ok]

    n = len(BATHY_LEVELS)
    lws = ([BATHY_LW * (1.0 + (BATHY_LW_DEEP_FACTOR - 1.0) * i / (n - 1)) for i in range(n)]
           if n > 1 else [BATHY_LW])
    cs = ax.contour(lon, lat, -elev, levels=BATHY_LEVELS, colors=BATHY_COLOR,
                    linewidths=lws, linestyles=BATHY_LS, transform=ccrs.PlateCarree(),
                    zorder=BATHY_ZORDER, alpha=BATHY_ALPHA)
    if BATHY_LABELS:
        ax.clabel(cs, levels=BATHY_LEVELS, fmt="%d m", fontsize=5, inline=True, inline_spacing=2)


def draw_eddies(ax, eddies):
    eddy_pe = ([pe.withStroke(linewidth=EDDY_HALO_LW, foreground=EDDY_HALO_COLOR,
                              alpha=EDDY_HALO_ALPHA)] if EDDY_HALO else None)
    for pts in eddies.get("ac_out", []) + eddies.get("ac_ms", []):
        ax.plot(*zip(*pts), "--", color=EDDY_COLOR, lw=EDDY_LW, zorder=EDDY_ZORDER,
                path_effects=eddy_pe, transform=ccrs.PlateCarree())
    for pts in eddies.get("cy_out", []) + eddies.get("cy_ms", []):
        ax.plot(*zip(*pts), "-", color=EDDY_COLOR, lw=EDDY_LW, zorder=EDDY_ZORDER,
                path_effects=eddy_pe, transform=ccrs.PlateCarree())
    for lo, la in eddies.get("ac_cen", []) + eddies.get("cy_cen", []):
        ax.plot(lo, la, "+", color=EDDY_COLOR, ms=EDDY_CENTER_MS, mew=1.4,
                zorder=EDDY_ZORDER + 1, path_effects=eddy_pe, transform=ccrs.PlateCarree())


def draw_currents(fig, ax, LON, LAT, u, v, is_zoom):
    """Geostrophic current vectors: uniform thin arrows, length = speed."""
    for nm, arr in (("ugos", u), ("vgos", v)):
        if np.ndim(arr) != 2:
            raise SystemExit(f"{nm} is {np.ndim(arr)}-D, expected 2-D.")
    s = QUIVER_STEP_ZOOM if is_zoom else QUIVER_STEP
    q_scale = QUIVER_SCALE_ZOOM if is_zoom else QUIVER_SCALE
    q_width = QUIVER_WIDTH_ZOOM if is_zoom else QUIVER_WIDTH
    MLON, MLAT = np.meshgrid(LON, LAT)
    q = ax.quiver(MLON[::s, ::s], MLAT[::s, ::s], u[::s, ::s], v[::s, ::s],
                  color=QUIVER_COLOR, alpha=QUIVER_ALPHA, transform=ccrs.PlateCarree(),
                  scale=q_scale, scale_units="width", width=q_width,
                  headwidth=QUIVER_HEADW, headlength=QUIVER_HEADL,
                  headaxislength=QUIVER_HEADL * 0.9, zorder=4)
    ax.quiverkey(q, 0.86, 1.02, QUIVER_KEY_SPEED, f"{QUIVER_KEY_SPEED:g} m s$^{{-1}}$",
                 labelpos="E", coordinates="axes", fontproperties={"size": 8})


def draw_ctds(ax, md_ctd, ag_ctd):
    for (lons, lats), col in ((md_ctd, MD_CTD_COLOR), (ag_ctd, AGU_CTD_COLOR)):
        if lons.size:
            ax.plot(lons, lats, "o", ls="None", color=col, ms=CTD_MS, mec=CTD_MEC,
                    mew=CTD_MEW, zorder=CTD_ZORDER, transform=ccrs.PlateCarree())


def setup_ticks(ax, extent, major=None):
    import matplotlib.ticker as mticker
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    lo0, lo1, la0, la1 = extent
    if major is None:
        major = TICK_MAJOR_ZOOM if (lo1 - lo0) <= 10 else TICK_MAJOR_FULL
    xt = np.arange(np.floor(lo0 / major) * major, np.ceil(lo1 / major) * major + 1e-6, major)
    yt = np.arange(np.floor(la0 / major) * major, np.ceil(la1 / major) * major + 1e-6, major)
    ax.set_xticks(xt, crs=ccrs.PlateCarree())
    ax.set_yticks(yt, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(N_MINOR))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(N_MINOR))
    ax.tick_params(which="major", length=6, width=0.8, direction="out",
                   top=True, right=True, labelsize=9)
    ax.tick_params(which="minor", length=3, width=0.5, direction="out", top=True, right=True)
    ax.tick_params(labeltop=False, labelright=False)
    ax.grid(which="major", color=GRID_COLOR, alpha=GRID_ALPHA, linewidth=GRID_LW,
            linestyle=":", zorder=6)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")


def make_figure(extent, suffix, date_stem, LON, LAT, adt, u, v, eddies, md_ctd, ag_ctd, cmap):
    """One map: ADT + 0.1 m isolines + currents + eddies + CTD positions."""
    figsize = FIG_SIZE_ZOOM if suffix else FIG_SIZE_FULL
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(list(extent), crs=ccrs.PlateCarree())

    if np.ndim(adt) != 2:
        raise SystemExit(f"adt is {np.ndim(adt)}-D, expected 2-D.")

    img = ax.imshow(adt, extent=[LON.min(), LON.max(), LAT.min(), LAT.max()],
                    cmap=cmap, origin="lower", vmin=ADT_VMIN, vmax=ADT_VMAX,
                    alpha=0.85, transform=ccrs.PlateCarree(), zorder=1)
    cb = fig.colorbar(img, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label("ADT (m)")

    draw_currents(fig, ax, LON, LAT, u, v, is_zoom=bool(suffix))
    draw_adt_isolines(ax, LON, LAT, adt)

    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, zorder=3)
    draw_bathymetry(ax, extent)
    draw_eddies(ax, eddies)
    draw_ctds(ax, md_ctd, ag_ctd)
    setup_ticks(ax, extent)

    handles = [
        mlines.Line2D([], [], color=ADT_ISO_COLOR, ls="-", lw=ADT_ISO_LW * 2,
                      label=f"ADT isolines ({ADT_ISO_STEP:g} m)"),
        mlines.Line2D([], [], color=EDDY_COLOR, ls="-", lw=EDDY_LW, label="Cyclone"),
        mlines.Line2D([], [], color=EDDY_COLOR, ls="--", lw=EDDY_LW, label="Anticyclone"),
        mlines.Line2D([], [], color=BATHY_COLOR, ls="-", lw=BATHY_LW * 2, label="Bathymetry"),
    ]
    if md_ctd[0].size:
        handles.append(mlines.Line2D([], [], color=MD_CTD_COLOR, ls="None", marker="o",
                                     ms=CTD_MS, mec=CTD_MEC, mew=CTD_MEW,
                                     label=f"MD CTD ({md_ctd[0].size})"))
    if ag_ctd[0].size:
        handles.append(mlines.Line2D([], [], color=AGU_CTD_COLOR, ls="None", marker="o",
                                     ms=CTD_MS, mec=CTD_MEC, mew=CTD_MEW,
                                     label=f"SA Agulhas CTD ({ag_ctd[0].size})"))
    ax.legend(handles=handles, loc="lower left", fontsize=7, framealpha=0.9)

    zoom_txt = "  -  zoom" if suffix else ""
    ax.set_title(f"WHIRLS  ADT + {ADT_ISO_STEP:g} m isolines + currents "
                 f"+ CTD stations  -  {date_stem}{zoom_txt}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{date_stem}_WHIRLS_adt_isolines_ctd{suffix}.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close(fig)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    target = argv[0] if argv else TARGET_DATE

    LON, LAT, adt, ugos, vgos, date_stem = ssh_for_date(target)
    eddies = eddies_for_date(date_stem)

    md_ctd = load_ctd_positions_nc(MD_CTD_NC, label="MD CTD")
    ag_ctd = load_ctd_positions_xlsx(find_agulhas_xlsx(), sheet=AGU_CTD_SHEET,
                                     dedupe_by_station=True, label="SA Agulhas CTD")

    cmap = getattr(cmc, ADT_CMAP_NAME)
    make_figure((LON_MIN, LON_MAX, LAT_MIN, LAT_MAX), "", date_stem,
                LON, LAT, adt, ugos, vgos, eddies, md_ctd, ag_ctd, cmap)
    make_figure(ZOOM_EXTENT, "_zoom", date_stem,
                LON, LAT, adt, ugos, vgos, eddies, md_ctd, ag_ctd, cmap)


if __name__ == "__main__":
    main()