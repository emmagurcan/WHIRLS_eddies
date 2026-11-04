"""
WHIRLS_adt_isolines_ctd.py
==========================

Static map of the ADT field for a chosen TARGET DATE (default 2026-07-13)
with its ADT ISOLINES every 0.1 m, the geostrophic current vectors, the
TOEddies contours, and the COMPLETE set of CTD station positions:

    Marion Dufresne (MD) CTDs   small ORANGE bullets
    SA Agulhas (AG) CTDs        small VIOLET bullets

The CTD positions are read from two Excel files (paths in the CONFIG block
below, variables MD_CTD_XLSX and AGU_CTD_XLSX):

    WHIRLS_MD_CTD_cast_metadata.xlsx
        sheet "CTD cast metadata"; one row per cast, columns
        "Lat [S]" / "Lon [E]" in signed decimal degrees.

    Agu073_WHIRLS_Cruise_Jun2026_CTD_profile_data_For_Distribution.xlsx
        sheet "Profile"; one row per DEPTH BIN, so the rows are
        de-duplicated by station to one position per CTD.

Both loaders match the lat/lon (and station) columns case-insensitively, so
minor header changes are tolerated. Reading .xlsx needs the openpyxl package
(pip install openpyxl).

The SSH file is the one in SSH/ whose date is NEAREST the target date (a
warning is printed if it is not the exact day). Eddies are taken for the
plotted SSH date, as in WHIRLS_adt_tracks_full_v2.py.

Figures produced:
    <date>_WHIRLS_adt_isolines_ctd.png          full region
    <date>_WHIRLS_adt_isolines_ctd_zoom.png     zoom

Usage (from the same working directory as the other WHIRLS routines, the
one that contains SSH/, TOEDDIES/, gebco_capetown.nc ...):

    python WHIRLS_adt_isolines_ctd.py                 target date 2026-07-13
    python WHIRLS_adt_isolines_ctd.py 2026-07-15      another target date

Data layout is otherwise the same as WHIRLS_adt_tracks_full_v2.py.
"""

import os
import sys
import glob
import json
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap

import cartopy.crs as ccrs
import cartopy.feature as cfeature


# ===========================================================================
#  CONFIG
# ===========================================================================
# Date of the field to plot: the SSH/<date>_cmems_data.nc NEAREST this date
# is used. Can be overridden from the command line (first argument).
TARGET_DATE = "2026-07-13"

LON_MIN, LON_MAX = 6.0, 22.0
LAT_MIN, LAT_MAX = -42.0, -30.0
ZOOM_EXTENT = (8.0, 16.0, -40.0, -34.0)

SSH_DIR  = "SSH"
EDDY_DIR = "TOEDDIES"

# ---------------------------------------------------------------------------
#  CTD position files  <-- EDIT THESE if the .xlsx files live elsewhere.
#
#  As written, both files are expected in the SAME directory you launch the
#  script from (the usual WHIRLS working directory). To keep them elsewhere,
#  give a relative or absolute path, e.g.:
#      MD_CTD_XLSX = os.path.join("CTD", "WHIRLS_MD_CTD_cast_metadata.xlsx")
#      MD_CTD_XLSX = "/Users/speich/WHIRLS/CTD/WHIRLS_MD_CTD_cast_metadata.xlsx"
# ---------------------------------------------------------------------------
MD_CTD_XLSX  = "/Users/speich/SABRINA/PROJETS/WHIRLS_SWOT/DATA/DATA_MD252/WHIRLS_MD_CTD_cast_metadata.xlsx"
##AGU_CTD_XLSX = ("/Users/speich/SABRINA/PROJETS/WHIRLS_SWOT/DATA/DATA_SA_Agulhas/Agu073_WHIRLS_Cruise_Jun2026_CTD_profile_data"
AGU_CTD_XLSX = "/Users/speich/SABRINA/PROJETS/WHIRLS_SWOT/DATA/DATA_SA_Agulhas/Agu073WHIRLSCruiseJun2026CTDprofiledata.xlsx"
MD_CTD_SHEET  = None      # None = first sheet ("CTD cast metadata")
AGU_CTD_SHEET = "Profile"

# CTD bullet styling: small filled circles, MD orange / Agulhas violet,
# with a thin dark edge so they read over any ADT colour.
MD_CTD_COLOR   = "#FF8C00"    # orange
AGU_CTD_COLOR  = "#9400D3"    # violet
CTD_MS         = 3.5          # small bullets
CTD_MEW        = 0.3
CTD_MEC        = "k"
CTD_ZORDER     = 8            # above eddies (5) and isolines (3)

GEBCO_FILE   = "gebco_capetown.nc"
BATHY_LEVELS = [100, 200, 500, 1000, 2000, 3000, 4000, 5000]
BATHY_COLOR  = "0.45"
BATHY_LW     = 0.35
BATHY_LW_DEEP_FACTOR = 1.25
BATHY_LS     = "solid"
BATHY_ALPHA  = 0.7
BATHY_ZORDER = 2.5
BATHY_LABELS = False

# --- Eddy contours: bold black, white halo; cyclones solid, anticyclones
#     dashed (same styling as WHIRLS_adt_isolines_tracks_full_v2.py) ---
EDDY_COLOR      = "k"
EDDY_LW         = 1.3
EDDY_ZORDER     = 5
EDDY_HALO       = True
EDDY_HALO_COLOR = "white"
EDDY_HALO_LW    = 2.4
EDDY_HALO_ALPHA = 0.85
EDDY_CENTER_MS  = 6

ADT_VMIN, ADT_VMAX = -0.7, 1.7
ADT_CMAP_NAME = "vik"
SCM_DIR = "/Users/speich/SABRINA2/COLORMAPS/ScientificColourMaps8"

# --- ADT isolines: contours of the plotted ADT field itself, every 0.1 m ---
ADT_ISO_STEP     = 0.1
ADT_ISO_COLOR    = "0.15"
ADT_ISO_LW       = 0.5
ADT_ISO_ALPHA    = 0.85
ADT_ISO_ZORDER   = 3.0
ADT_ISO_LABELS   = True
ADT_ISO_LABEL_FS = 5

QUIVER_STEP       = 2
QUIVER_STEP_ZOOM  = 1
QUIVER_SCALE      = 110
QUIVER_SCALE_ZOOM = 45
QUIVER_WIDTH      = 0.0007
QUIVER_WIDTH_ZOOM = 0.0012
QUIVER_HEADW      = 3.0
QUIVER_HEADL      = 3.5
QUIVER_ALPHA      = 0.9
QUIVER_COLOR      = "k"
QUIVER_KEY_SPEED  = 1.0

# --- Axes ticks, subticks, grid, and output resolution ---
TICK_MAJOR_FULL = 2.0
TICK_MAJOR_ZOOM = 1.0
N_MINOR         = 4
GRID_COLOR      = "white"
GRID_ALPHA      = 0.45
GRID_LW         = 0.4
FIG_DPI         = 300
FIG_SIZE_FULL   = (13, 9)
FIG_SIZE_ZOOM   = (12, 9)
# ===========================================================================


def adt_colormap():
    name = ADT_CMAP_NAME
    for p in (os.path.join(SCM_DIR, name, f"{name}.txt"),
              os.path.join(SCM_DIR, name, f"{name}.dat"),
              os.path.join(SCM_DIR, f"{name}.txt")):
        if os.path.exists(p):
            rgb = np.loadtxt(p)
            if rgb.max() > 1.0:
                rgb = rgb / 255.0
            print(f"[info] ADT colormap '{name}' from {p}")
            return LinearSegmentedColormap.from_list(name, rgb)
    try:
        from cmcrameri import cm as cmc
        return getattr(cmc, name)
    except Exception:
        print(f"[info] Crameri '{name}' not found; using RdBu_r.")
        return plt.get_cmap("RdBu_r")


def load_geojson(path):
    if not os.path.exists(path):
        return {"features": []}
    with open(path) as f:
        return json.load(f)


def polygon_coords(gj):
    out = []
    for feat in gj.get("features", []):
        g = feat.get("geometry")
        if not g or not g.get("coordinates"):
            continue
        if g["type"] == "Polygon":
            out.append(g["coordinates"][0])
        elif g["type"] == "MultiPolygon":
            out.append(g["coordinates"][0][0])
    return out


def point_coords(gj):
    return [tuple(f["geometry"]["coordinates"][:2])
            for f in gj.get("features", [])
            if f.get("geometry") and f["geometry"].get("coordinates")]


def _subset_bounds(ds, name, lo, hi):
    """Slice a coordinate to [lo, hi] whatever its storage order."""
    vals = ds[name].values
    if vals.size > 1 and vals[0] > vals[-1]:      # descending
        return {name: slice(hi, lo)}
    return {name: slice(lo, hi)}


def ssh_for_date(target):
    """
    Open the SSH file whose date is NEAREST `target` (a 'YYYY-MM-DD' string);
    return (LON, LAT, adt, ugos, vgos, date_stem). Warns when the nearest
    available file is not the exact target day.
    """
    files = sorted(glob.glob(os.path.join(SSH_DIR, "*_cmems_data.nc")))
    if not files:
        raise SystemExit(f"No SSH files in {SSH_DIR}/")

    try:
        tgt = datetime.strptime(target, "%Y-%m-%d")
    except ValueError:
        raise SystemExit(f"TARGET_DATE '{target}' is not YYYY-MM-DD.")

    def _stem(p):
        return os.path.basename(p).split("_cmems_data")[0]

    def _gap(p):
        try:
            return abs((datetime.strptime(_stem(p), "%Y-%m-%d") - tgt).days)
        except ValueError:
            return 10**6

    path = min(files, key=_gap)
    stem = _stem(path)
    gap = _gap(path)
    if gap == 0:
        print(f"SSH file for {target}: {path}")
    else:
        print(f"[warn] no SSH file for {target}; using the nearest, "
              f"{stem} ({gap} day{'s' if gap != 1 else ''} away): {path}")

    ds = xr.open_dataset(path)
    sel = {}
    sel.update(_subset_bounds(ds, "longitude", LON_MIN, LON_MAX))
    sel.update(_subset_bounds(ds, "latitude", LAT_MIN, LAT_MAX))
    sub = ds.sel(**sel)

    nlat, nlon = sub.sizes.get("latitude", 0), sub.sizes.get("longitude", 0)
    if nlat < 2 or nlon < 2:
        full_lat = ds["latitude"].values
        full_lon = ds["longitude"].values
        ds.close()
        raise SystemExit(
            f"\nSSH subset is empty or degenerate: {nlat} lat x {nlon} lon.\n"
            f"  requested window : lon {LON_MIN}..{LON_MAX}, "
            f"lat {LAT_MIN}..{LAT_MAX}\n"
            f"  file covers      : lon {full_lon.min():.3f}.."
            f"{full_lon.max():.3f}, lat {full_lat.min():.3f}.."
            f"{full_lat.max():.3f}\n"
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

    print(f"  grid {LAT.size} lat x {LON.size} lon; "
          f"lon {LON.min():.3f}..{LON.max():.3f}, "
          f"lat {LAT.min():.3f}..{LAT.max():.3f}")
    return LON, LAT, adt, u, v, stem


def eddies_for_date(ssh_date):
    """Eddy geometry for the plotted SSH date (nearest file, with warning)."""
    cens = sorted(glob.glob(os.path.join(
        EDDY_DIR, "*_anticyclones_centers.geojson")))
    if not cens:
        print(f"[warn] no eddy files in {EDDY_DIR}/; eddies skipped.")
        return {}

    available = [os.path.basename(p).split("_anticyclones")[0] for p in cens]
    if ssh_date in available:
        date = ssh_date
        print(f"Eddies from {date} (matches the SSH date)")
    else:
        def _days(d):
            try:
                return abs((datetime.strptime(d, "%Y-%m-%d") -
                            datetime.strptime(ssh_date, "%Y-%m-%d")).days)
            except ValueError:
                return 10**6
        date = min(available, key=_days)
        gap = _days(date)
        print(f"[warn] no eddy file for the SSH date {ssh_date}; using "
              f"{date} instead ({gap} day{'s' if gap != 1 else ''} away).")

    def gj(kind, prod):
        return load_geojson(os.path.join(
            EDDY_DIR, f"{date}_{kind}_{prod}.geojson"))
    return {
        "ac_out": polygon_coords(gj("anticyclones", "outcontour")),
        "cy_out": polygon_coords(gj("cyclones", "outcontour")),
        "ac_ms":  polygon_coords(gj("anticyclones", "maxspeed")),
        "cy_ms":  polygon_coords(gj("cyclones", "maxspeed")),
        "ac_cen": point_coords(gj("anticyclones", "centers")),
        "cy_cen": point_coords(gj("cyclones", "centers")),
    }


# ---------------------------------------------------------------------------
#  CTD positions from the two Excel files
# ---------------------------------------------------------------------------
def _excel_serial_to_date(x):
    """Excel serial day number -> datetime (None if not parseable)."""
    try:
        return datetime(1899, 12, 30) + timedelta(days=float(x))
    except (TypeError, ValueError):
        return None


def load_ctd_positions_xlsx(path, sheet=None, dedupe_by_station=False,
                            label="CTD"):
    """
    Read CTD positions from an .xlsx file; return (lons, lats) as arrays.

    The header row is found automatically (first row containing a cell whose
    text starts with 'lat'), and the latitude / longitude / station / date
    columns are matched case-insensitively ('Lat [S]', 'latitude', ... all
    work). With `dedupe_by_station`, consecutive rows of the same station
    (one row per depth bin, as in the Agulhas Profile sheet) collapse to a
    single position per station.

    If the header says [S] (or [W]) but every value is positive, the sign is
    flipped, with a warning.
    """
    if not os.path.exists(path):
        print(f"[warn] {label}: file not found: {path}; no positions drawn. "
              f"Edit MD_CTD_XLSX / AGU_CTD_XLSX at the top of this script.")
        return np.array([]), np.array([])
    try:
        from openpyxl import load_workbook
    except ImportError:
        raise SystemExit(
            "Reading .xlsx needs the openpyxl package. Install it with:\n"
            "    python -m pip install openpyxl")

    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet] if (sheet and sheet in wb.sheetnames) else wb[wb.sheetnames[0]]
    if sheet and sheet not in wb.sheetnames:
        print(f"[warn] {label}: sheet '{sheet}' not in {path}; using "
              f"'{ws.title}'. Sheets: {', '.join(wb.sheetnames)}")

    rows = ws.iter_rows(values_only=True)

    # locate the header row and the columns
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
        raise SystemExit(f"{label}: no latitude column found in {path}.")

    def find(*prefixes):
        for i, c in enumerate(header):
            if any(c..startswith(p) for p in prefixes):
                return i
        return None

    i_lon = find("lon", "long")
    i_stn = find("station", "cast", "point")
    i_date = find("date", "mm/dd")
    if i_lon is None:
        wb.close()
        raise SystemExit(f"{label}: no longitude column found in {path}.")

    lats, lons, dates, last_key = [], [], [], object()
    for r in rows:
        if r is None or max(i_lat, i_lon) >= len(r):
            continue
        try:
            la = float(r[i_lat])
            lo = float(r[i_lon])
        except (TypeError, ValueError):
            continue
        key = (r[i_stn] if (i_stn is not None and i_stn < len(r)) else
               (round(la, 4), round(lo, 4)))
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

    lats = np.array(lats, dtype=float)
    lons = np.array(lons, dtype=float)

    # hemisphere fix: '[s]'/'[w]' headers with all-positive values
    if lats.size and "[s]" in header[i_lat] and (lats > 0).all():
        print(f"[warn] {label}: latitude header says [S] but values are "
              f"positive; flipping sign.")
        lats = -lats
    if lons.size and "[w]" in header[i_lon] and (lons > 0).all():
        print(f"[warn] {label}: longitude header says [W] but values are "
              f"positive; flipping sign.")
        lons = -lons

    span = ""
    if dates:
        span = (f"; dates {min(dates):%Y-%m-%d} .. {max(dates):%Y-%m-%d}")
    print(f"{label}: {lats.size} CTD positions from {path}"
          f" (sheet '{ws.title}'){span}")
    return lons, lats


# ---------------------------------------------------------------------------
#  Drawing
# ---------------------------------------------------------------------------
def draw_adt_isolines(ax, LON, LAT, adt):
    """Contour the plotted ADT field every ADT_ISO_STEP metres (levels snap
    to whole multiples of the step over the field's own range)."""
    finite = np.isfinite(adt)
    if not finite.any():
        print("[warn] ADT field is all-NaN; isolines skipped.")
        return None
    lo = np.floor(np.nanmin(adt) / ADT_ISO_STEP) * ADT_ISO_STEP
    hi = np.ceil(np.nanmax(adt) / ADT_ISO_STEP) * ADT_ISO_STEP
    levels = np.arange(lo, hi + ADT_ISO_STEP / 2.0, ADT_ISO_STEP)
    cs = ax.contour(LON, LAT, adt, levels=levels, colors=ADT_ISO_COLOR,
                    linewidths=ADT_ISO_LW, linestyles="solid",
                    transform=ccrs.PlateCarree(), zorder=ADT_ISO_ZORDER,
                    alpha=ADT_ISO_ALPHA)
    if ADT_ISO_LABELS:
        ax.clabel(cs, fmt="%.1f", fontsize=ADT_ISO_LABEL_FS,
                  inline=True, inline_spacing=2)
    print(f"ADT isolines: {levels.size} levels, "
          f"{levels[0]:.1f} .. {levels[-1]:.1f} m, step {ADT_ISO_STEP:g} m")
    return cs


def draw_bathymetry(ax, extent):
    if not os.path.exists(GEBCO_FILE):
        print(f"[warn] bathymetry not found: {GEBCO_FILE}")
        return
    lo0, lo1, la0, la1 = extent
    ds = xr.open_dataset(GEBCO_FILE)
    dv = next((v for v in ("elevation", "z", "Band1", "depth")
               if v in ds.variables), None)
    if dv is None:
        ds.close()
        return
    latn = "lat" if "lat" in ds.coords else "latitude"
    lonn = "lon" if "lon" in ds.coords else "longitude"
    lat, lon = ds[latn].values, ds[lonn].values
    elev = np.squeeze(ds[dv].values).astype(float)
    ds.close()
    la_ok = (lat >= la0 - 1) & (lat <= la1 + 1)
    lo_ok = (lon >= lo0 - 1) & (lon <= lo1 + 1)
    lat, lon, elev = lat[la_ok], lon[lo_ok], elev[np.ix_(la_ok, lo_ok)]

    n = len(BATHY_LEVELS)
    if n > 1:
        lws = [BATHY_LW * (1.0 + (BATHY_LW_DEEP_FACTOR - 1.0) * i / (n - 1))
               for i in range(n)]
    else:
        lws = [BATHY_LW]

    cs = ax.contour(lon, lat, -elev, levels=BATHY_LEVELS, colors=BATHY_COLOR,
                    linewidths=lws, linestyles=BATHY_LS,
                    transform=ccrs.PlateCarree(), zorder=BATHY_ZORDER,
                    alpha=BATHY_ALPHA)
    if BATHY_LABELS:
        ax.clabel(cs, levels=BATHY_LEVELS, fmt="%d m", fontsize=5,
                  inline=True, inline_spacing=2)


def draw_eddies(ax, eddies):
    import matplotlib.patheffects as pe
    eddy_pe = ([pe.withStroke(linewidth=EDDY_HALO_LW,
                              foreground=EDDY_HALO_COLOR,
                              alpha=EDDY_HALO_ALPHA)] if EDDY_HALO else None)
    for pts in eddies.get("ac_out", []) + eddies.get("ac_ms", []):
        ax.plot(*zip(*pts), "--", color=EDDY_COLOR, lw=EDDY_LW,
                zorder=EDDY_ZORDER, path_effects=eddy_pe,
                transform=ccrs.PlateCarree())
    for pts in eddies.get("cy_out", []) + eddies.get("cy_ms", []):
        ax.plot(*zip(*pts), "-", color=EDDY_COLOR, lw=EDDY_LW,
                zorder=EDDY_ZORDER, path_effects=eddy_pe,
                transform=ccrs.PlateCarree())
    for lo, la in eddies.get("ac_cen", []) + eddies.get("cy_cen", []):
        ax.plot(lo, la, "+", color=EDDY_COLOR, ms=EDDY_CENTER_MS, mew=1.4,
                zorder=EDDY_ZORDER + 1, path_effects=eddy_pe,
                transform=ccrs.PlateCarree())


def draw_currents(fig, ax, LON, LAT, u, v, is_zoom):
    """Geostrophic current vectors: uniform thin arrows, length = speed."""
    for nm, arr in (("ugos", u), ("vgos", v)):
        if np.ndim(arr) != 2:
            raise SystemExit(
                f"{nm} is {np.ndim(arr)}-D, expected 2-D. The velocity field "
                f"was lost between reading and plotting. Value: {arr!r}")

    s = QUIVER_STEP_ZOOM if is_zoom else QUIVER_STEP
    q_scale = QUIVER_SCALE_ZOOM if is_zoom else QUIVER_SCALE
    q_width = QUIVER_WIDTH_ZOOM if is_zoom else QUIVER_WIDTH

    MLON, MLAT = np.meshgrid(LON, LAT)
    q = ax.quiver(MLON[::s, ::s], MLAT[::s, ::s], u[::s, ::s], v[::s, ::s],
                  color=QUIVER_COLOR, alpha=QUIVER_ALPHA,
                  transform=ccrs.PlateCarree(), scale=q_scale,
                  scale_units="width", width=q_width,
                  headwidth=QUIVER_HEADW, headlength=QUIVER_HEADL,
                  headaxislength=QUIVER_HEADL * 0.9, zorder=4)
    ax.quiverkey(q, 0.86, 1.02, QUIVER_KEY_SPEED,
                 f"{QUIVER_KEY_SPEED:g} m s$^{{-1}}$", labelpos="E",
                 coordinates="axes", fontproperties={"size": 8})


def draw_ctds(ax, md_ctd, ag_ctd):
    """Small bullets at every CTD position: MD orange, Agulhas violet."""
    if md_ctd[0].size:
        ax.plot(md_ctd[0], md_ctd[1], "o", ls="None",
                color=MD_CTD_COLOR, ms=CTD_MS, mec=CTD_MEC, mew=CTD_MEW,
                zorder=CTD_ZORDER, transform=ccrs.PlateCarree())
    if ag_ctd[0].size:
        ax.plot(ag_ctd[0], ag_ctd[1], "o", ls="None",
                color=AGU_CTD_COLOR, ms=CTD_MS, mec=CTD_MEC, mew=CTD_MEW,
                zorder=CTD_ZORDER, transform=ccrs.PlateCarree())


def setup_ticks(ax, extent, major=None):
    import matplotlib.ticker as mticker
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    lo0, lo1, la0, la1 = extent
    if major is None:
        major = TICK_MAJOR_ZOOM if (lo1 - lo0) <= 10 else TICK_MAJOR_FULL

    xt = np.arange(np.floor(lo0 / major) * major,
                   np.ceil(lo1 / major) * major + 1e-6, major)
    yt = np.arange(np.floor(la0 / major) * major,
                   np.ceil(la1 / major) * major + 1e-6, major)
    ax.set_xticks(xt, crs=ccrs.PlateCarree())
    ax.set_yticks(yt, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(N_MINOR))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(N_MINOR))
    ax.tick_params(which="major", length=6, width=0.8, direction="out",
                   top=True, right=True, labelsize=9)
    ax.tick_params(which="minor", length=3, width=0.5, direction="out",
                   top=True, right=True)
    ax.tick_params(labeltop=False, labelright=False)
    ax.grid(which="major", color=GRID_COLOR, alpha=GRID_ALPHA,
            linewidth=GRID_LW, linestyle=":", zorder=6)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")


def make_figure(extent, suffix, date_stem, LON, LAT, adt, u, v,
                eddies, md_ctd, ag_ctd):
    """One map: ADT + 0.1 m isolines + currents + eddies + CTD positions."""
    cmap = adt_colormap()
    figsize = FIG_SIZE_ZOOM if suffix else FIG_SIZE_FULL
    fig, ax = plt.subplots(figsize=figsize,
                           subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(list(extent), crs=ccrs.PlateCarree())

    if np.ndim(adt) != 2:
        raise SystemExit(
            f"adt is {np.ndim(adt)}-D, expected 2-D. The SSH field was lost "
            f"between reading and plotting. Value: {adt!r}")

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
        mlines.Line2D([], [], color=ADT_ISO_COLOR, ls="-",
                      lw=ADT_ISO_LW * 2,
                      label=f"ADT isolines ({ADT_ISO_STEP:g} m)"),
        mlines.Line2D([], [], color=EDDY_COLOR, ls="-", lw=EDDY_LW,
                      label="Cyclone"),
        mlines.Line2D([], [], color=EDDY_COLOR, ls="--", lw=EDDY_LW,
                      label="Anticyclone"),
        mlines.Line2D([], [], color=BATHY_COLOR, ls="-", lw=BATHY_LW * 2,
                      label="Bathymetry"),
    ]
    if md_ctd[0].size:
        handles.append(mlines.Line2D(
            [], [], color=MD_CTD_COLOR, ls="None", marker="o", ms=CTD_MS,
            mec=CTD_MEC, mew=CTD_MEW,
            label=f"MD CTD ({md_ctd[0].size})"))
    if ag_ctd[0].size:
        handles.append(mlines.Line2D(
            [], [], color=AGU_CTD_COLOR, ls="None", marker="o", ms=CTD_MS,
            mec=CTD_MEC, mew=CTD_MEW,
            label=f"SA Agulhas CTD ({ag_ctd[0].size})"))
    ax.legend(handles=handles, loc="lower left", fontsize=7, framealpha=0.9)

    zoom_txt = "  -  zoom" if suffix else ""
    ax.set_title(f"WHIRLS  ADT + {ADT_ISO_STEP:g} m isolines + currents "
                 f"+ CTD stations  -  {date_stem}{zoom_txt}")

    out = f"{date_stem}_WHIRLS_adt_isolines_ctd{suffix}.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close(fig)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    target = argv[0] if argv else TARGET_DATE

    LON, LAT, adt, ugos, vgos, date_stem = ssh_for_date(target)
    eddies = eddies_for_date(date_stem)

    md_ctd = load_ctd_positions_xlsx(
        MD_CTD_XLSX, sheet=MD_CTD_SHEET, dedupe_by_station=False,
        label="MD CTD")
    ag_ctd = load_ctd_positions_xlsx(
        AGU_CTD_XLSX, sheet=AGU_CTD_SHEET, dedupe_by_station=True,
        label="SA Agulhas CTD")

    make_figure((LON_MIN, LON_MAX, LAT_MIN, LAT_MAX), "", date_stem,
                LON, LAT, adt, ugos, vgos, eddies, md_ctd, ag_ctd)
    make_figure(ZOOM_EXTENT, "_zoom", date_stem,
                LON, LAT, adt, ugos, vgos, eddies, md_ctd, ag_ctd)


if __name__ == "__main__":
    main()
