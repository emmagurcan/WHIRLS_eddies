#!/usr/bin/env python
"""Track a WHIRLS eddy through the TOEddies segment graph and plot its trajectory.

Merging / splitting events come from the TOEddies graph (not decided by hand). The eddy
is followed backwards as far as the graph goes (and forwards up to --end), keeping the
lowest-cost branch at each event. Eddies that split off are drawn in gray; eddies that
merged in can be drawn too (--show-merged).

  mamba activate toeddies_pkgs
  python analysis/plot_eddy_lineage.py AE11 --dataset AVISO
  python analysis/plot_eddy_lineage.py CE18 --dataset VARDYN
  python analysis/plot_eddy_lineage.py all  --dataset all           # 2 eddies x 2 datasets
  # any other eddy:
  python analysis/plot_eddy_lineage.py MYEDDY --type AC --seed 2026-07-05 11.3 -37.5

Outputs in figures/eddy_tracks/toeddies_lineage/:
  lineage_<NAME>_<DATASET>.png           the map
  lineage_<NAME>_<DATASET>_track.csv     main track (one row per detection)
  lineage_<NAME>_<DATASET>_events.csv    merging / splitting events
  lineage_<NAME>_<DATASET>_branches.csv  split-off (and merged-in) eddies
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
from matplotlib.colors import Normalize
from matplotlib.path import Path as MplPath
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, str(Path(__file__).resolve().parent))
import eddy_lineage as el  # noqa: E402

BASE = el.BASE
OUT_DIR = BASE / "figures/eddy_tracks/toeddies_lineage"
FNAME_BATHY = BASE / "data/bathymetry/etopo2.nc"
SADCP_TRACK_CSV = BASE / "figures/sadcp_vs_geostrophy/sadcp_vs_geostrophy_track.csv"
DIR_LADCP = BASE / "data/LADCP/traitement_avec_SADCP"

# The WHIRLS eddies. The seed is a day when the ship sampled the eddy; its position is
# taken from the manual track (csv) on that day, then matched to a TOEddies eddy.
EDDIES = {
    "AE11": dict(name="AE11_20260514", orientation="Anticyclones", seed_date="2026-07-05",
                 csv=BASE / "figures/eddy_tracks/tracked_eddy_AE11_20260514.csv"),
    "CE18": dict(name="CE18_20260601", orientation="Cyclones", seed_date="2026-07-16",
                 csv=BASE / "figures/eddy_tracks/tracked_eddy_CE18_202601601.csv"),
}

EVENT_STYLE = {"merging": dict(marker="P", color="magenta"),
               "splitting": dict(marker="X", color="cyan")}
# light gray reads well on the dark-blue bathymetry and is distinct from the red/yellow track
BRANCH_STYLE = {"split": dict(color="0.88", ls="-", lw=1.3),
                "merge": dict(color="0.7", ls="--", lw=1.1)}
crs = ccrs.PlateCarree()
ft = 10


# =============================================================================
def seed_from_csv(csv, date):
    df = pd.read_csv(csv, parse_dates=["date"])
    row = df.iloc[(df.date - pd.Timestamp(date)).abs().argmin()]
    return row.date + pd.Timedelta(hours=12), float(row.lon), float(row.lat)


def subsample(df, every_days):
    """Indices of the rows to draw a contour for: one every `every_days` days."""
    if not every_days or df.empty:
        return []
    t = (df.time - df.time.iloc[0]) / pd.Timedelta(days=every_days)
    k = np.floor(t + 1e-6).astype(int)
    idx = list(df.index[~k.duplicated()])
    if idx[-1] != df.index[-1]:
        idx.append(df.index[-1])
    return idx


def setup_map(fig, extent):
    proj = ccrs.Mercator(central_longitude=(extent[0] + extent[1]) / 2)
    ax = fig.add_subplot(111, projection=proj)
    ax.set_extent(extent, crs=crs)
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="gray", zorder=3)
    ax.coastlines(resolution="10m", color="k", linewidth=0.6, zorder=4)
    gl = ax.gridlines(linewidth=1, color="black", alpha=0.2, linestyle="--",
                      draw_labels=True, zorder=6)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": ft}
    x0, x1, y0, y1 = extent
    topo = xr.open_dataset(FNAME_BATHY)["topo"].sel(lon=slice(x0 - 1, x1 + 1),
                                                     lat=slice(y0 - 1, y1 + 1))
    c = ax.contourf(topo.lon, topo.lat, topo, np.arange(-5500, 1, 500),
                    transform=crs, cmap="Blues_r", zorder=0, alpha=0.85, extend="min")
    return ax, c


def in_situ(main):
    """SADCP positions / LADCP stations inside the eddy's outermost contour that day."""
    day_xy = {}
    for _, r in main.iterrows():
        xy = r.outer if r.outer is not None else r.contour
        if xy is not None:
            day_xy.setdefault(r.time.normalize(), xy)
    inside = lambda xy, x, y: MplPath(xy).contains_points(np.c_[x, y])
    ship = pd.DataFrame(columns=["time", "lon", "lat"])
    try:
        df = pd.read_csv(SADCP_TRACK_CSV, parse_dates=["time"]).dropna(subset=["lon", "lat"])
        keep = [g[inside(day_xy[d], g.lon, g.lat)]
                for d, g in df.groupby(df.time.dt.normalize()) if d in day_xy]
        if keep:
            ship = pd.concat(keep)
    except Exception as e:
        print(f"  SADCP not loaded ({e})")
    st_in = pd.DataFrame(columns=["station", "lon", "lat"])
    try:
        from ladcp_tools import load_all_stations
        st = load_all_stations(str(DIR_LADCP) + "/")
        tcol = "datetime" if "datetime" in st.columns else st.select_dtypes("datetime").columns[0]
        ok = [i for i, s in st.iterrows()
              if pd.Timestamp(s[tcol]).normalize() in day_xy
              and inside(day_xy[pd.Timestamp(s[tcol]).normalize()], [s.lon], [s.lat])[0]]
        st_in = st.loc[ok]
    except Exception as e:
        print(f"  LADCP not loaded ({e})")
    return ship, st_in


def plot_lineage(lin, name, extent=None, pad=1.0, contour_every=3, center_every=1,
                 show_merged=False, insitu=False, out=None):
    main = lin.main.reset_index(drop=True)
    typ = "Anticyclone" if lin.seed["orientation"] == "Anticyclones" else "Cyclone"

    if extent is None:   # everything that is drawn: main track, branches, their contours
        drawn = [main] + lin.splits + (lin.merges if show_merged else [])
        cont = [c for d in drawn for c in d.contour if c is not None]
        xs = np.concatenate([d.lon.values for d in drawn] + [c[:, 0] for c in cont])
        ys = np.concatenate([d.lat.values for d in drawn] + [c[:, 1] for c in cont])
        extent = [xs.min() - pad, xs.max() + pad, ys.min() - pad, ys.max() + pad]

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax, cbathy = setup_map(fig, extent)
    tnum = mdates.date2num(main.time)
    norm = Normalize(vmin=tnum.min(), vmax=max(tnum.max(), tnum.min() + 1))
    cmap = plt.get_cmap("autumn")

    # --- other eddies (gray) --------------------------------------------------------
    handles = []
    kinds = [("split", lin.splits)] + ([("merge", lin.merges)] if show_merged else [])
    for kind, branches in kinds:
        st = BRANCH_STYLE[kind]
        for b in branches:
            b = b.reset_index(drop=True)
            for k in subsample(b, contour_every):
                xy = b.contour[k]
                if xy is not None:
                    ax.plot(xy[:, 0], xy[:, 1], color=st["color"], ls=st["ls"],
                            lw=st["lw"] * 0.9, alpha=0.8, transform=crs, zorder=7)
            ax.plot(b.lon, b.lat, color=st["color"], ls=st["ls"], lw=st["lw"],
                    transform=crs, zorder=7.5)
            bc = b.loc[subsample(b, center_every)]
            ax.scatter(bc.lon, bc.lat, s=8, c=st["color"], edgecolors="k", linewidths=0.2,
                       transform=crs, zorder=7.5)
            end = b.iloc[-1] if kind == "split" else b.iloc[0]
            ax.annotate(f"{end.time:%d %b}", (end.lon, end.lat), xycoords=crs._as_mpl_transform(ax),
                        xytext=(3, -8), textcoords="offset points", fontsize=5.5,
                        color=st["color"], zorder=7.5)
        if branches:
            lbl = ("Eddy that split off (track, contours)" if kind == "split"
                   else "Eddy that merged in (track, contours)")
            handles.append(mlines.Line2D([], [], marker=".", **st, label=lbl))

    # --- main track ---------------------------------------------------------------
    ax.plot(main.lon, main.lat, "-", color="k", lw=0.8, transform=crs, zorder=9)
    for k in subsample(main, contour_every):
        xy = main.contour[k]
        if xy is not None:
            ax.plot(xy[:, 0], xy[:, 1], color=cmap(norm(tnum[k])), lw=1.2,
                    transform=crs, zorder=8)
    daily = main.loc[subsample(main, center_every)]   # one dot every center_every days
    ax.scatter(daily.lon, daily.lat, c=mdates.date2num(daily.time), cmap=cmap, norm=norm,
               s=14, edgecolors="k", linewidths=0.3, transform=crs, zorder=10)
    s = lin.seed
    first = main.iloc[0]
    ax.scatter(first.lon, first.lat, marker="o", s=70, facecolors="none", edgecolors="k",
               linewidths=1.5, transform=crs, zorder=11)
    handles = [
        mlines.Line2D([], [], marker="o", ls="-", color="k", mfc="orange", ms=5,
                      label="Center" + (" (daily)" if center_every == 1
                                        else f" (every {center_every:g} d)")),
        mlines.Line2D([], [], color="orange", lw=1.2,
                      label=f"Max-speed contour (every {contour_every:g} d)"),
        mlines.Line2D([], [], marker="o", ls="", mfc="none", mec="k", mew=1.5, ms=8,
                      label=f"First detection ({first.time:%d %b}"
                            + (", start of data)" if first.time <= s["data_start"] else ")")),
    ] + handles

    # --- in-situ ------------------------------------------------------------------
    if insitu:
        ship, st = in_situ(main)
        if len(ship):
            ax.scatter(ship.lon, ship.lat, s=2, c="k", transform=crs, zorder=12)
            handles.append(mlines.Line2D([], [], marker=".", ls="", color="k",
                                         label="Ship (SADCP) inside eddy"))
        if len(st):
            ax.scatter(st.lon, st.lat, marker="^", s=40, c="lime", edgecolors="k",
                       transform=crs, zorder=13)
            handles.append(mlines.Line2D([], [], marker="^", ls="", mfc="lime", mec="k",
                                         label="LADCP station inside eddy"))

    # --- events -------------------------------------------------------------------
    seen = set()
    for _, e in lin.events.iterrows():
        st = EVENT_STYLE[e.kind]
        ax.scatter(e.lon, e.lat, marker=st["marker"], s=130, c=st["color"], edgecolors="k",
                   linewidths=0.8, transform=crs, zorder=15)
        ax.annotate(f"{e.time:%d %b}", xy=(e.lon, e.lat), xycoords=crs._as_mpl_transform(ax),
                    xytext=(7, 7) if e.kind == "splitting" else (7, -14),
                    textcoords="offset points", fontsize=7, fontweight="bold",
                    zorder=16, bbox=dict(boxstyle="round,pad=0.2", fc="white",
                                         ec=st["color"], alpha=0.9))
        if e.kind not in seen:
            seen.add(e.kind)
            handles.append(mlines.Line2D([], [], marker=st["marker"], ls="", color=st["color"],
                                         mec="k", ms=9, label=f"{e.kind.capitalize()} event"))

    ax.legend(handles=handles, loc="best", fontsize=7, framealpha=0.9).set_zorder(20)
    ax.set_title(f"{typ} {name} ({s['dataset']}, TOEddies): "
                 f"{main.time.min():%b %d} – {main.time.max():%b %d, %Y}", fontsize=12)
    cax = inset_axes(ax, width="80%", height="3%", loc="lower center",
                     bbox_to_anchor=(0, -0.12, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
    fig.colorbar(cbathy, cax=cax, orientation="horizontal").set_label("Depth [m]", fontsize=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb2 = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.6, pad=0.03)
    cb2.ax.yaxis.set_major_locator(mdates.AutoDateLocator(maxticks=6))
    cb2.ax.yaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    cb2.set_label("Date", fontsize=9)
    if out:
        fig.savefig(out, bbox_inches="tight")
        print(f"  saved {out}")
    plt.close(fig)


def save_tables(lin, stem):
    cols = ["time", "lon", "lat", "r_maxspeed_km", "r_outer_km", "speed", "ssh",
            "run", "segment", "node"]
    lin.main[cols].to_csv(f"{stem}_track.csv", index=False)
    lin.events.drop(columns="contour").to_csv(f"{stem}_events.csv", index=False)
    br = [b[cols + ["role", "branch"]] for b in lin.splits + lin.merges]
    (pd.concat(br) if br else pd.DataFrame(columns=cols + ["role", "branch"])).to_csv(
        f"{stem}_branches.csv", index=False)


# =============================================================================
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("eddy", help="AE11, CE18, all, or a name for a custom --seed")
    p.add_argument("--dataset", default="AVISO", choices=["AVISO", "VARDYN", "all"])
    p.add_argument("--seed", nargs=3, metavar=("DATE", "LON", "LAT"),
                   help="seed for a custom eddy (needs --type)")
    p.add_argument("--seed-date", help="override the preset seed date (position from csv)")
    p.add_argument("--type", choices=["AC", "C"], help="for a custom eddy")
    p.add_argument("--start", default=None, help="do not go back before this date")
    p.add_argument("--end", default="2026-07-27",
                   help="stop the forward tracking at this date ('none' = end of life)")
    p.add_argument("--branch-days", type=float, default=15,
                   help="follow split / merged eddies at most this many days (0 = whole life)")
    p.add_argument("--contour-every", type=float, default=3, help="days between contours")
    p.add_argument("--show-merged", action="store_true",
                   help="also draw the eddies that merged in (dashed gray)")
    p.add_argument("--no-stitch", action="store_true",
                   help="VarDyn: do not stitch the may/june/july runs")
    p.add_argument("--center-every", type=float, default=1, help="days between centers")
    p.add_argument("--insitu", action="store_true", help="add SADCP / LADCP inside the eddy")
    p.add_argument("--extent", nargs=4, type=float, metavar=("W", "E", "S", "N"))
    a = p.parse_args()

    end = None if str(a.end).lower() == "none" else a.end
    if a.seed:
        if not a.type:
            p.error("--seed needs --type AC|C")
        todo = {a.eddy: dict(name=a.eddy, seed=(pd.Timestamp(a.seed[0]), float(a.seed[1]),
                                                float(a.seed[2])),
                             orientation="Anticyclones" if a.type == "AC" else "Cyclones")}
    else:
        keys = list(EDDIES) if a.eddy == "all" else [a.eddy]
        todo = {}
        for k in keys:
            cfg = dict(EDDIES[k])
            cfg["seed"] = seed_from_csv(cfg["csv"], a.seed_date or cfg["seed_date"])
            todo[k] = cfg
    datasets = list(el.DATASETS) if a.dataset == "all" else [a.dataset]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for k, cfg in todo.items():
        for ds in datasets:
            t, lon, lat = cfg["seed"]
            print(f"\n=== {cfg['name']} / {ds}: seed {t:%Y-%m-%d} ({lon:.2f}, {lat:.2f})")
            try:
                lin = el.track(ds, cfg["orientation"], t, lon, lat, t_min=a.start,
                               t_max=end, branch_days=a.branch_days or None, stitch=not a.no_stitch)
            except (ValueError, RuntimeError) as e:
                print("  FAILED:", e)
                continue
            m = lin.main
            print(f"  main track {m.time.min():%Y-%m-%d} -> {m.time.max():%Y-%m-%d}, "
                  f"{len(lin.events)} events, {len(lin.splits)} split-off, "
                  f"{len(lin.merges)} merged-in eddies")
            for _, e in lin.events.iterrows():
                print(f"    {e.time:%Y-%m-%d %Hh}  {e.kind:9s} ({e.lon:6.2f}, {e.lat:6.2f}) "
                      f"{e.n_branches} branches [{e.run}]")
            stem = OUT_DIR / f"lineage_{cfg['name']}_{ds}"
            save_tables(lin, stem)
            plot_lineage(lin, cfg["name"], extent=a.extent, contour_every=a.contour_every,
                         show_merged=a.show_merged, insitu=a.insitu,
                         center_every=a.center_every,
                         out=f"{stem}.png")


if __name__ == "__main__":
    matplotlib.use("Agg")
    main()
