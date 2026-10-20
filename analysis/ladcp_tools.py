"""
Tools to read LADCP .lad files, plot a single profile, and map all stations.

Usage (from a terminal):
    python ladcp_tools.py profile path/to/002.lad
    python ladcp_tools.py map path/to/traitement_avec_SADCP

Or import the functions in a notebook:
    from ladcp_tools import read_lad, plot_profile, load_all_stations, plot_station_map
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. Reading a .lad file
# ---------------------------------------------------------------------------

def _parse_coord(text):
    """Convert a string like "34°S 30.0003'" or "16°E 38.3016'" to decimal degrees."""
    m = re.search(r"(\d+)\D*?([NSEW])\s*([\d.]+)", text)
    if m is None:
        raise ValueError(f"Cannot parse coordinate: {text!r}")
    deg, hemi, minutes = int(m.group(1)), m.group(2), float(m.group(3))
    value = deg + minutes / 60.0
    return -value if hemi in "SW" else value


def read_lad(path):
    """
    Read one .lad file.

    Returns
    -------
    meta : dict   station info (station, date, time, lat, lon, deviation, ...)
    df   : DataFrame with columns z [m], u [m/s], v [m/s], ev [m/s], speed [m/s]
    """
    path = Path(path)
    # The degree sign may be stored in latin-1 or utf-8, so decode leniently
    lines = path.read_bytes().decode("latin-1").splitlines()

    meta = {"station": path.stem, "file": str(path)}
    columns = ["z", "u", "v", "ev"]
    data_start = 0

    for i, line in enumerate(lines):
        if "=" not in line:
            data_start = i
            break
        key, value = [s.strip() for s in line.split("=", 1)]
        if key == "Start_Lat":
            meta["lat"] = _parse_coord(value)
        elif key == "Start_Lon":
            meta["lon"] = _parse_coord(value)
        elif key == "Columns":
            columns = value.split(":")
        elif key == "Deviation":
            meta["deviation"] = float(value)   # magnetic declination used [deg]
        else:
            meta[key.lower()] = value

    rows = [l.split() for l in lines[data_start:] if l.strip()]
    df = pd.DataFrame(rows, columns=columns).astype(float)
    df["speed"] = np.hypot(df["u"], df["v"])

    meta["datetime"] = pd.to_datetime(
        f"{meta.get('date', '')} {meta.get('start_time', '')}", errors="coerce"
    )
    meta["max_depth"] = df["z"].max()
    return meta, df


# ---------------------------------------------------------------------------
# 2. Plotting one profile
# ---------------------------------------------------------------------------

def plot_profile(path, save=None):
    """Plot u, v (with error band) and speed versus depth for one .lad file."""
    meta, df = read_lad(path)
    z = df["z"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6), sharey=True, dpi=150)

    for comp, color, label in [("u", "tab:blue", "u (east)"), ("v", "tab:red", "v (north)")]:
        ax1.fill_betweenx(z, df[comp] - df["ev"], df[comp] + df["ev"],
                          color=color, alpha=0.2, linewidth=0)
        ax1.plot(df[comp], z, color=color, label=label)
    ax1.axvline(0, color="k", linewidth=0.6)
    ax1.set_xlabel("Velocity [m/s]")
    ax1.set_ylabel("Depth [m]")
    ax1.legend(loc="lower right", fontsize="small")
    ax1.set_title("Velocity components (±ev)", fontsize="medium")

    ax2.plot(df["speed"], z, color="k")
    ax2.set_xlabel("Speed [m/s]")
    ax2.set_xlim(left=0)
    ax2.set_title("Speed", fontsize="medium")

    ax1.invert_yaxis()
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)

    fig.suptitle(
        f"LADCP station {meta['station']}  |  {meta['datetime']:%Y-%m-%d %H:%M}  |  "
        f"{meta['lat']:.3f}°, {meta['lon']:.3f}°"
    )
    fig.tight_layout()
    if save:
        fig.savefig(save, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 3. Collecting all stations
# ---------------------------------------------------------------------------

def load_all_stations(root, layer=(0, 200)):
    """
    Read every .lad file under `root` (e.g. traitement_avec_SADCP/002/002.lad).

    `layer` is a depth range [m] used to compute a layer-mean velocity
    (u_mean, v_mean), handy to draw an arrow per station on the map.

    Returns a DataFrame with one row per station.
    """
    records = []
    for f in sorted(Path(root).glob("*/*.lad")):
        try:
            meta, df = read_lad(f)
        except Exception as e:
            print(f"Skipping {f}: {e}")
            continue
        sel = df[(df["z"] >= layer[0]) & (df["z"] <= layer[1])]
        meta["u_mean"] = sel["u"].mean()
        meta["v_mean"] = sel["v"].mean()
        meta["mean_ev"] = df["ev"].mean()
        records.append(meta)

    stations = pd.DataFrame(records).sort_values("station").reset_index(drop=True)
    print(f"Found {len(stations)} stations in {root}")
    return stations


# ---------------------------------------------------------------------------
# 4. Map of all stations
# ---------------------------------------------------------------------------

def plot_station_map(stations, extent=None, ax=None, crs=None, arrows=True,
                     arrow_scale=3, ft=6, fill_land=False, save=None):
    """
    Plot all LADCP stations on a map.

    If you pass an existing cartopy `ax` (e.g. your ADT map), stations are added
    on top of it. Otherwise a new map is created.

    Station marker colour = cast date. Arrows = layer-mean velocity
    (computed in load_all_stations).
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    crs = crs or ccrs.PlateCarree()

    if ax is None:
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(), zorder=3)
        # fill_land=True needs the Natural Earth land file; cartopy downloads it
        # the first time, which can hang on a slow or blocked network.
        if fill_land:
            ax.add_feature(cfeature.LAND.with_scale("50m"), linewidth=0.5, zorder=4,
                           edgecolor="black", facecolor="gray")
        ax.coastlines(resolution="10m", color="k", zorder=5)
        grid = ax.gridlines(linewidth=1, color="black", alpha=0.2,
                            linestyle="--", draw_labels=True, zorder=6)
        grid.top_labels = False
        grid.right_labels = False
        grid.xlabel_style = {"size": ft}
        grid.ylabel_style = {"size": ft}
    else:
        fig = ax.figure

    if extent is None:
        pad = 0.5
        extent = [stations.lon.min() - pad, stations.lon.max() + pad,
                  stations.lat.min() - pad, stations.lat.max() + pad]
    ax.set_extent(extent, crs=crs)

    # Colour stations by date so you can follow the cruise in time
    dates = stations["datetime"].map(pd.Timestamp.toordinal)
    sc = ax.scatter(stations.lon, stations.lat, c=dates, cmap="viridis",
                    s=15, edgecolors="black", linewidths=0.4,
                    zorder=100, transform=crs)

    for _, st in stations.iterrows():
        ax.annotate(
            str(int(st.station)) if st.station.isdigit() else st.station,
            xy=(st.lon, st.lat), xycoords=crs._as_mpl_transform(ax),
            xytext=(3, 3), textcoords="offset points",
            fontsize=4, ha="left", va="bottom", zorder=101,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                      edgecolor="black", linewidth=0.3, alpha=0.8),
        )

    if arrows and "u_mean" in stations:
        Q = ax.quiver(stations.lon, stations.lat, stations.u_mean, stations.v_mean,
                      transform=crs, scale=arrow_scale, width=0.003, zorder=99)
        ax.quiverkey(Q, X=0.85, Y=1.03, U=0.5, label="0.5 m/s",
                     labelpos="E", fontproperties={"size": ft})

    cbar = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    ticks = np.linspace(dates.min(), dates.max(), 5).round()
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([pd.Timestamp.fromordinal(int(t)).strftime("%d %b") for t in ticks])
    cbar.ax.tick_params(labelsize=ft)
    cbar.set_label("Cast date", fontsize=ft)

    if save:
        fig.savefig(save, bbox_inches="tight")
    return fig, ax

def load_eddies(dir_eddies):
    """Load cyclonic and anticyclonic eddy shapefiles."""

    anti = {}
    cyc = {}

    for fname in os.listdir(dir_eddies):

        if "anti" in fname:
            target = anti
        else:
            target = cyc

        if "maxspeed" in fname:
            target["maxspeed"] = gpd.read_file(os.path.join(dir_eddies, fname))

        elif "outcontour" in fname:
            target["outcontour"] = gpd.read_file(os.path.join(dir_eddies, fname))

        elif "centers" in fname:
            target["centers"] = gpd.read_file(os.path.join(dir_eddies, fname))

    return anti, cyc
# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("profile", "map"):
        print(__doc__)
        sys.exit(1)

    if sys.argv[1] == "profile":
        f = Path(sys.argv[2])
        plot_profile(f, save=f"ladcp_profile_{f.stem}.png")
        print(f"Saved ladcp_profile_{f.stem}.png")
    else:
        st = load_all_stations(sys.argv[2])
        st.drop(columns=["file"]).to_csv("ladcp_stations.csv", index=False)
        plot_station_map(st, save="ladcp_stations_map.png")
        print("Saved ladcp_stations_map.png and ladcp_stations.csv")