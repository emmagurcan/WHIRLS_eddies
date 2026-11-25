"""
WHIRLS: for every CTD cast, check whether it was done inside a TOEddies eddy contour
(on the day of the cast) or not. The casts outside every eddy are the candidates for
a cruise-based background T / S / density profile.

Input  : data/CTD/ctd_metadata.nc (cast positions + time_start), TOEddies shapefiles
Output : data/CTD/ctd_toeddies_classification.csv (one row per cast)
         figures/ctd_toeddies_classification.png  (map of the casts)

Put this file next to tools_eg.py (it uses get_eddy_files / load_eddies from there).
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from shapely.geometry import Point

from tools_eg import BASE, get_eddy_files, load_eddies

# ===========================================================================
#  CONFIG
# ===========================================================================
DIR_CTD = BASE / "data/CTD"
MD_CTD_NC = DIR_CTD / "ctd_metadata.nc"
OUT_CSV = DIR_CTD / "ctd_toeddies_classification.csv"
OUT_FIG = BASE / "figures" / "ctd_toeddies_classification.png"

CONTOUR = "outcontour"     # "outcontour" (outer edge of the eddy) or "maxspeed"
BUFFER_KM = 0.0            # > 0: casts closer than this to a contour also count as "in eddy"
R_KM = 6371.0
# ===========================================================================


# ---------------------------------------------------------------------------
#  CTD metadata (same variable search as the cruise-report notebook)
# ---------------------------------------------------------------------------
def _find_position_var(ds, short, std_name, units):
    """Name of the lat (or lon) variable of the cast START position."""
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
        return (("end" in nl) or ("bottom" in nl), nl not in (short, std_name), "start" not in nl)
    return sorted(cands, key=rank)[0]


def load_ctd_metadata(path):
    """-> DataFrame(cast, station, time, lon, lat), one row per cast."""
    ds = xr.open_dataset(path)
    print(f"{path.name}: variables {list(ds.variables)}")
    lat_n = _find_position_var(ds, "lat", "latitude", "degrees_north")
    lon_n = _find_position_var(ds, "lon", "longitude", "degrees_east")
    if lat_n is None or lon_n is None:
        raise SystemExit("No latitude / longitude variable found in the CTD metadata.")
    time_n = "time_start" if "time_start" in ds.variables else next(
        (n for n, v in ds.variables.items() if np.issubdtype(v.dtype, np.datetime64)), None)
    if time_n is None:
        raise SystemExit("No time variable found in the CTD metadata (expected 'time_start').")
    stn_n = next((n for n in ds.variables
                  if any(k in n.lower() for k in ("station", "stn", "cast", "profile"))
                  and n not in (lat_n, lon_n, time_n)), None)
    print(f"  using lat='{lat_n}', lon='{lon_n}', time='{time_n}', station='{stn_n}'")

    lat = np.asarray(ds[lat_n].values, float).ravel()
    lon = np.asarray(ds[lon_n].values, float).ravel()
    time = pd.to_datetime(np.asarray(ds[time_n].values).ravel())
    if stn_n is not None:
        stn = np.asarray(ds[stn_n].values).ravel()
        stn = [s.decode().strip() if isinstance(s, bytes) else s for s in stn]
    else:
        stn = list(range(len(lat)))
    ds.close()
    if not (len(lat) == len(lon) == len(time) == len(stn)):
        raise SystemExit(f"lat/lon/time/station have different lengths: "
                         f"{len(lat)}, {len(lon)}, {len(time)}, {len(stn)}")
    df = pd.DataFrame(dict(cast=np.arange(len(lat)), station=stn, time=time, lon=lon, lat=lat))
    bad = df[["lon", "lat"]].isna().any(axis=1) | df["time"].isna()
    if bad.any():
        print(f"  ! {bad.sum()} cast(s) without position or time skipped")
    return df[~bad].reset_index(drop=True)


# ---------------------------------------------------------------------------
#  TOEddies contours of one day (cached)
# ---------------------------------------------------------------------------
_eddy_cache = {}


def contours_for_day(day):
    """day 'YYYY-MM-DD' -> (list of (type, index, polygon), product name) for CONTOUR."""
    if day in _eddy_cache:
        return _eddy_cache[day]
    files, source = get_eddy_files(day.replace("-", ""), day)
    polys = []
    if files:
        anti, cyc = load_eddies(files)
        for kind, eddies in (("AC", anti), ("C", cyc)):
            gdf = eddies.get(CONTOUR)
            if gdf is None:
                continue
            for i, geom in enumerate(gdf.geometry):
                if geom is not None and not geom.is_empty:
                    polys.append((kind, i, geom))
    _eddy_cache[day] = (polys, source)
    return polys, source


def dist_to_contour_km(lon, lat, geom):
    """Distance (km) from the cast to the nearest vertex of the contour."""
    ring = geom.exterior if geom.geom_type == "Polygon" else max(geom.geoms, key=lambda g: g.area).exterior
    xy = np.asarray(ring.coords)
    dlat = np.deg2rad(xy[:, 1] - lat)
    dlon = np.deg2rad(xy[:, 0] - lon)
    a = np.sin(dlat / 2) ** 2 + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(xy[:, 1])) * np.sin(dlon / 2) ** 2
    return float(np.min(2 * R_KM * np.arcsin(np.sqrt(a))))


# ---------------------------------------------------------------------------
#  Classification
# ---------------------------------------------------------------------------
def classify(ctd):
    rows = []
    for _, c in ctd.iterrows():
        day = c["time"].strftime("%Y-%m-%d")
        polys, source = contours_for_day(day)
        row = dict(c, eddy_product=source, status="no TOEddies file", eddy_type="",
                   eddy_index=-1, nearest_type="", dist_nearest_km=np.nan)
        if source is None:
            rows.append(row)
            continue
        pt = Point(c["lon"], c["lat"])
        best = (np.inf, "", -1)
        inside = None
        for kind, i, geom in polys:
            if geom.covers(pt):
                inside = (kind, i)
                best = (0.0, kind, i)
                break
            d = dist_to_contour_km(c["lon"], c["lat"], geom)
            if d < best[0]:
                best = (d, kind, i)
        row.update(dist_nearest_km=best[0] if np.isfinite(best[0]) else np.nan,
                   nearest_type=best[1])
        if inside is not None:
            row.update(status=f"inside {inside[0]}", eddy_type=inside[0], eddy_index=inside[1])
        elif best[0] <= BUFFER_KM:
            row.update(status=f"within {BUFFER_KM:g} km of {best[1]}", eddy_type=best[1],
                       eddy_index=best[2])
        else:
            row["status"] = "outside"
        rows.append(row)
    out = pd.DataFrame(rows)
    out["outside_eddies"] = out["status"] == "outside"
    return out


def plot_casts(df):
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="lightgray", zorder=2)
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6, zorder=3)
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.4, linestyle="--")
        gl.top_labels = gl.right_labels = False
        tr = dict(transform=ccrs.PlateCarree())
    except ImportError:
        fig, ax = plt.subplots(figsize=(9, 7))
        tr = {}
    style = {"outside": ("o", "tab:green"), "inside AC": ("^", "tab:red"),
             "inside C": ("v", "tab:blue"), "no TOEddies file": ("x", "0.5")}
    handles = []
    for status, g in df.groupby("status"):
        m, col = style.get(status, ("s", "orange"))
        ax.plot(g["lon"], g["lat"], m, ls="None", color=col, ms=6, mec="k", mew=0.4, zorder=5, **tr)
        handles.append(mlines.Line2D([], [], marker=m, ls="None", color=col, mec="k", mew=0.4,
                                     label=f"{status} ({len(g)})"))
    pad = 0.5
    lims = [df["lon"].min() - pad, df["lon"].max() + pad, df["lat"].min() - pad, df["lat"].max() + pad]
    if tr:
        ax.set_extent(lims, crs=tr["transform"])
    else:
        ax.set_xlim(lims[:2])
        ax.set_ylim(lims[2:])
    ax.legend(handles=handles, loc="lower left", fontsize=8, framealpha=0.9)
    ax.set_title(f"WHIRLS CTD casts vs TOEddies {CONTOUR} (day of each cast)", fontsize=11)
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    print(f"Saved {OUT_FIG}")
    return fig


def main():
    ctd = load_ctd_metadata(MD_CTD_NC)
    print(f"{len(ctd)} casts, {ctd['time'].min():%Y-%m-%d} .. {ctd['time'].max():%Y-%m-%d}")
    df = classify(ctd)

    print("\n" + df["status"].value_counts().to_string())
    print("\nCasts outside every eddy contour:")
    cols = ["cast", "station", "time", "lon", "lat", "dist_nearest_km", "nearest_type", "eddy_product"]
    print(df.loc[df["outside_eddies"], cols].to_string(index=False,
          formatters={"lon": "{:.3f}".format, "lat": "{:.3f}".format,
                      "dist_nearest_km": "{:.0f}".format}))

    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {OUT_CSV}")
    plot_casts(df)
    plt.show()
    return df


if __name__ == "__main__":
    df = main()
