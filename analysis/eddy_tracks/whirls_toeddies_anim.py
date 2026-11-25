"""
WHIRLS: daily animation of the TOEddies contours on ADT, 11-21 July 2026.

Same map as the cruise-report figure (cruise_report_figure.ipynb, cell 7, panel background):
ADT colour fill (Crameri 'vik'), geostrophic current arrows, ADT isolines every 0.1 m,
land + coastline, bathymetry contours, TOEddies (cyclones solid, anticyclones dashed,
black with a white halo, '+' at the centres). No instruments.

One PNG per day goes to FRAME_DIR, then ffmpeg joins them into OUT_MP4
(or a GIF if ffmpeg is not found). Days without an SSH file are skipped.
"""
import re
import shutil
import subprocess
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ===========================================================================
#  CONFIG
# ===========================================================================
START_DATE = "2026-06-27"
END_DATE   = "2026-07-27"

BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
DIR_SSH = BASE / "data/SSH"
DIR_EDDIES = {"TOEDDIES_better_ADT": BASE / "data/TOEDDIES/TOEDDIES_better_ADT",
              "TOEDDIES": BASE / "data/TOEDDIES/TOEDDIES"}
EDDY_PRODUCT = "TOEDDIES_better_ADT"        # tried first; the other one is the fallback
FNAME_BATHY = BASE / "data/bathymetry/etopo2.nc"

OUT_DIR = BASE / "figures/toeddies_animation"
FRAME_DIR = OUT_DIR / "frames"
OUT_MP4 = OUT_DIR / f"WHIRLS_TOEddies_ADT_{START_DATE}_{END_DATE}.mp4"
FPS = 2                                     # frames (days) per second
ANIM_DPI = 150                              # per frame (the report figure uses 300)

# Map area (same as the cruise-report figure)
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = 7.5, 15.0, -40.0, -33.0
extent = (LON_MIN, LON_MAX, LAT_MIN, LAT_MAX)

SHOW_CURRENTS = True                        # geostrophic arrows, as in the report figure

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

ADT_VMIN, ADT_VMAX = -0.5, 1.5
ADT_CMAP_NAME = "vik"       # Crameri 'vik' if cmcrameri is installed, else cmocean balance
ADT_ALPHA = 0.8

# ADT isolines every 0.1 m
ADT_ISO_STEP = 0.1
ADT_ISO_COLOR = "0.15"
ADT_ISO_LW = 0.5
ADT_ISO_ALPHA = 0.85
ADT_ISO_ZORDER = 3.0
ADT_ISO_LABELS = True
ADT_ISO_LABEL_FS = 5

QUIVER_STEP = 2
QUIVER_SCALE = 110
QUIVER_WIDTH = 0.0007
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

FIG_SIZE = (9, 9)           # fixed, so every frame has the same size
LEG_KW = dict(loc="upper center", bbox_to_anchor=(0.5, -0.07), fontsize=7, frameon=False,
              columnspacing=1.4, handlelength=2.4, ncol=4)
TITLE_KW = dict(loc="left", fontsize=11, pad=16)
# ===========================================================================


def adt_colormap():
    try:
        from cmcrameri import cm as cmc
        print(f"[info] ADT colormap: Crameri '{ADT_CMAP_NAME}'")
        return getattr(cmc, ADT_CMAP_NAME)
    except Exception:
        pass
    try:
        import cmocean
        print("[info] cmcrameri not installed; ADT colormap: cmocean 'balance'")
        return cmocean.cm.balance
    except Exception:
        print("[info] ADT colormap: RdBu_r")
        return plt.get_cmap("RdBu_r")


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
    """SSH file in DIR_SSH for exactly `target` (YYYY-MM-DD) ->
    (LON, LAT, adt, ugos, vgos, 'YYYY-MM-DD'), or None if there is no file that day."""
    tgt = datetime.strptime(target, "%Y-%m-%d")
    files = [f for f in sorted(DIR_SSH.glob("*.nc")) if file_date(f) == tgt]
    if not files:
        print(f"[warn] no SSH file for {target}; day skipped")
        return None
    path = files[0]
    print(f"SSH file for {target}: {path}")

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
    return LON, LAT, adt, u, v, target


# ---------------------------------------------------------------------------
#  TOEddies (shapefiles or geojson, read with geopandas)
# ---------------------------------------------------------------------------
def _geom_coords(geom):
    """Polygon -> exterior ring; MultiPolygon -> first polygon; Point -> (x, y)."""
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
    """Eddy geometry for the SSH date (nearest available date, with warning).
    Returns (eddies dict, label 'product, YYYY-MM-DD')."""
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
        name = "TOEddies (better ADT)" if product == "TOEDDIES_better_ADT" else "TOEddies"
        return out, f"{name}, {best:%Y-%m-%d}"

    print("[warn] no TOEddies files found; eddies skipped.")
    return {}, None


# ---------------------------------------------------------------------------
#  Bathymetry (read once)
# ---------------------------------------------------------------------------
def load_bathymetry(extent):
    if not FNAME_BATHY.exists():
        print(f"[warn] bathymetry not found: {FNAME_BATHY}")
        return None
    lo0, lo1, la0, la1 = extent
    ds = xr.open_dataset(FNAME_BATHY)
    dv = next((v for v in ("topo", "elevation", "z", "Band1", "depth") if v in ds.variables), None)
    if dv is None:
        ds.close()
        print(f"[warn] no elevation variable in {FNAME_BATHY}; bathymetry skipped.")
        return None
    latn = "lat" if "lat" in ds.coords else "latitude"
    lonn = "lon" if "lon" in ds.coords else "longitude"
    lat, lon = ds[latn].values, ds[lonn].values
    la_ok = (lat >= la0 - 1) & (lat <= la1 + 1)
    lo_ok = (lon >= lo0 - 1) & (lon <= lo1 + 1)
    elev = ds[dv].isel({latn: np.flatnonzero(la_ok), lonn: np.flatnonzero(lo_ok)})
    elev = np.squeeze(elev.transpose(latn, lonn).values).astype(float)
    ds.close()
    return lon[lo_ok], lat[la_ok], elev


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
    return cs


def draw_bathymetry(ax, bathy):
    if bathy is None:
        return
    lon, lat, elev = bathy
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


def draw_currents(ax, LON, LAT, u, v):
    """Geostrophic current vectors: uniform thin arrows, length = speed."""
    s = QUIVER_STEP
    MLON, MLAT = np.meshgrid(LON, LAT)
    q = ax.quiver(MLON[::s, ::s], MLAT[::s, ::s], u[::s, ::s], v[::s, ::s],
                  color=QUIVER_COLOR, alpha=QUIVER_ALPHA, transform=ccrs.PlateCarree(),
                  scale=QUIVER_SCALE, scale_units="width", width=QUIVER_WIDTH,
                  headwidth=QUIVER_HEADW, headlength=QUIVER_HEADL,
                  headaxislength=QUIVER_HEADL * 0.9, zorder=4)
    ax.quiverkey(q, 0.86, 1.04, QUIVER_KEY_SPEED, f"{QUIVER_KEY_SPEED:g} m s$^{{-1}}$",
                 labelpos="E", coordinates="axes", fontproperties={"size": 8})


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


def background_handles():
    return [
        mlines.Line2D([], [], color=ADT_ISO_COLOR, ls="-", lw=ADT_ISO_LW * 2, alpha=ADT_ISO_ALPHA,
                      label=f"ADT isolines ({ADT_ISO_STEP:g} m)"),
        mlines.Line2D([], [], color=EDDY_COLOR, ls="-", lw=EDDY_LW, label="Cyclone"),
        mlines.Line2D([], [], color=EDDY_COLOR, ls="--", lw=EDDY_LW, label="Anticyclone"),
        mlines.Line2D([], [], color=BATHY_COLOR, ls="-", lw=BATHY_LW * 2, label="Bathymetry"),
    ]


# ---------------------------------------------------------------------------
#  One frame
# ---------------------------------------------------------------------------
def make_frame(day, cmap, bathy):
    ssh = ssh_for_date(day)
    if ssh is None:
        return None
    LON, LAT, adt, ugos, vgos, date_stem = ssh
    eddies, eddy_label = eddies_for_date(date_stem)

    fig = plt.figure(figsize=FIG_SIZE)
    ax = fig.add_axes([0.10, 0.16, 0.74, 0.76], projection=ccrs.PlateCarree())
    ax.set_extent(list(extent), crs=ccrs.PlateCarree())

    img = ax.imshow(adt, extent=[LON.min(), LON.max(), LAT.min(), LAT.max()],
                    cmap=cmap, origin="lower", vmin=ADT_VMIN, vmax=ADT_VMAX,
                    alpha=ADT_ALPHA, transform=ccrs.PlateCarree(), zorder=1,
                    interpolation="bilinear")
    if SHOW_CURRENTS:
        draw_currents(ax, LON, LAT, ugos, vgos)
    draw_adt_isolines(ax, LON, LAT, adt)
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="lightgray", zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6, zorder=3)
    draw_bathymetry(ax, bathy)
    draw_eddies(ax, eddies)

    setup_ticks(ax, extent)
    ax.set_extent(list(extent), crs=ccrs.PlateCarree())   # ticks must not widen the map

    ax.legend(handles=background_handles(), **LEG_KW)
    ax.set_title(f"ADT and TOEddies: {datetime.strptime(date_stem, '%Y-%m-%d'):%-d %B %Y}",
                 **TITLE_KW)
    if eddy_label:
        ax.text(0.99, 0.01, eddy_label, transform=ax.transAxes, fontsize=6, ha="right",
                va="bottom", zorder=200,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    # colour bar on the right, placed by hand from the map's real size (same on every frame)
    fig.canvas.draw()
    p = ax.get_position()
    cax = fig.add_axes([p.x1 + 0.015, p.y0 + 0.12 * p.height, 0.015, 0.76 * p.height])
    cb = fig.colorbar(img, cax=cax)
    cb.set_label("ADT (m)")
    cb.solids.set_alpha(1)

    out = FRAME_DIR / f"{date_stem}_TOEddies_ADT.png"
    fig.savefig(out, dpi=ANIM_DPI)          # no bbox_inches="tight": identical frame sizes
    plt.close(fig)
    print(f"Saved {out}")
    return out


def make_movie(frames):
    if shutil.which("ffmpeg"):
        lst = FRAME_DIR / "frames.txt"
        lst.write_text("".join(f"file '{f.resolve()}'\nduration {1 / FPS}\n" for f in frames)
                       + f"file '{frames[-1].resolve()}'\n")
        cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(lst),
               "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2,fps=25",
               "-c:v", "libx264", "-pix_fmt", "yuv420p", "-preset", "veryfast", str(OUT_MP4)]
        subprocess.run(cmd, check=True)
        print(f"Saved {OUT_MP4}")
    else:
        from PIL import Image
        gif = OUT_MP4.with_suffix(".gif")
        ims = [Image.open(f).convert("RGB") for f in frames]
        ims[0].save(gif, save_all=True, append_images=ims[1:], duration=int(1000 / FPS), loop=0)
        print(f"ffmpeg not found; saved {gif}")


def main():
    FRAME_DIR.mkdir(parents=True, exist_ok=True)
    cmap = adt_colormap()
    bathy = load_bathymetry(extent)

    d0 = datetime.strptime(START_DATE, "%Y-%m-%d")
    d1 = datetime.strptime(END_DATE, "%Y-%m-%d")
    days = [(d0 + timedelta(days=k)).strftime("%Y-%m-%d") for k in range((d1 - d0).days + 1)]

    frames = [f for f in (make_frame(day, cmap, bathy) for day in days) if f is not None]
    print(f"{len(frames)} frames ({len(days) - len(frames)} days skipped)")
    if frames:
        make_movie(frames)


if __name__ == "__main__":
    main()