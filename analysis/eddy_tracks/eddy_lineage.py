"""Follow one eddy through the TOEddies segment graph (merging / splitting included).

TOEddies (Trajectories/*.py) produces, for each orientation, a graph of *segments*
(``graph_AC/output_final.gt``, ``graph_C/output_final.gt``):

  * a vertex = a segment = a chain of instantaneous eddies with no event inside
    (``vp.inst_eddies`` lists their node ids, node = date_index * e_overestim + eddy_index);
  * an edge  = an overlap between the last eddy of one segment and the first eddy of the
    next one, with a ``cost`` (lower = more similar eddies);
  * in-degree  >= 2  -> MERGING  at the first eddy of the segment;
  * out-degree >= 2  -> SPLITTING at the last eddy of the segment.

Starting from a seed (time, lon, lat), we find the instantaneous eddy there, its segment,
and walk the graph backwards (and forwards). At every event we keep the branch with the
lowest cost as "the" eddy (same choice as TOEddies' own trajectories.py) and remember the
other branches: eddies that split off (followed forward) and eddies that merged in
(followed backward).

Several TOEddies runs covering successive, overlapping periods (VarDyn may/june/july) are
stitched: when the track reaches the overlap with the previous / next run, it is re-seeded
in that run at the middle of the overlap.

Usage: see plot_eddy_lineage.py.
"""

import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
sys.path.insert(0, str(BASE / "TOEddies/Detection_eddies/Trajectories"))
import util_eddies  # noqa: E402  (TOEddies' own helpers)

R_EARTH_KM = 6371.0
GRAPH_DIR = {"Anticyclones": "graph_AC", "Cyclones": "graph_C"}

# Date index of TOEddies -> time.
#   AVISO : days since 1950-01-01 (CNES julian day)
#   VarDyn: 6-h steps since 2026-01-01 00h (see run_toeddies_vardyn.sh)
DATASETS = {
    "AVISO": dict(dirs=["AVISO_BIG"], epoch="1950-01-01", step_h=24),
    "VARDYN": dict(dirs=["VARDYN_may", "VARDYN_june", "VARDYN_july"],
                   epoch="2026-01-01", step_h=6),
}
DIR_TOEDDIES = BASE / "data/TOEDDIES"


def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, (lon1, lat1, lon2, lat2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(a))


def _points_in_poly(xy, x, y):
    from matplotlib.path import Path as MplPath
    return MplPath(xy).contains_point((x, y))


# =============================================================================
# One TOEddies run (one directory, one orientation)
# =============================================================================
class Run:
    def __init__(self, run_dir, orientation, epoch, step_h):
        import graph_tool

        self.dir = Path(run_dir)
        self.name = self.dir.name
        self.orientation = orientation
        self.epoch = pd.Timestamp(epoch)
        self.step = pd.Timedelta(hours=step_h)

        gfile = self.dir / GRAPH_DIR[orientation] / "output_final.gt"
        if not gfile.exists():
            raise FileNotFoundError(f"{gfile} not found (TOEddies not finished?)")
        self.g = graph_tool.load_graph(str(gfile))
        if self.g.gp.orientation != orientation:
            raise FileNotFoundError(f"{gfile} is a graph of {self.g.gp.orientation}, "
                                    f"not {orientation} (check cyclone= in eddy_graph and the "
                                    f"orientation argument of graph_properties.py)")
        self.e_over = int(self.g.gp.e_overestim)
        self.shpc = util_eddies.SHPC_class(str(self.dir))
        self.d_first = self.shpc.d_init[0]
        self.d_last = self.shpc.d_init[-1] + int(self.shpc.n_dates[-1]) - 1
        self.t_first, self.t_last = self.time(self.d_first), self.time(self.d_last)

        self.node2seg = {}
        for v in self.g.iter_vertices():
            for n in self.g.vp.inst_eddies[v]:
                self.node2seg[int(n)] = int(v)

        self.cost = {(int(s), int(t)): float(c)
                     for s, t, c in self.g.get_edges(eprops=[self.g.ep.cost])}

    def __repr__(self):
        return (f"Run({self.name}, {self.orientation}, "
                f"{self.t_first:%Y-%m-%d %H}h -> {self.t_last:%Y-%m-%d %H}h)")

    # --- time ------------------------------------------------------------------
    def time(self, date_index):
        return self.epoch + int(date_index) * self.step

    def date_index(self, t):
        return int(round((pd.Timestamp(t) - self.epoch) / self.step))

    # --- graph -----------------------------------------------------------------
    def nodes(self, seg):
        return [int(n) for n in self.g.vp.inst_eddies[seg]]

    def preds(self, seg):
        """[(pred_segment, cost)] sorted by cost."""
        out = [(int(u), self.cost[(int(u), seg)]) for u in self.g.get_in_neighbors(seg)]
        return sorted(out, key=lambda x: x[1])

    def succs(self, seg):
        out = [(int(w), self.cost[(seg, int(w))]) for w in self.g.get_out_neighbors(seg)]
        return sorted(out, key=lambda x: x[1])

    # --- instantaneous eddies ----------------------------------------------------
    @lru_cache(maxsize=None)
    def eddy(self, node):
        """Center, max-speed contour (outermost if no max-speed contour) and outermost
        contour of one instantaneous eddy."""
        d, i = util_eddies.node_to_date_eddy(node, self.e_over)
        i_slice, ishape = self.shpc.comp_ishape(d, i, self.orientation)
        ext = self.shpc.get_reader(i_slice, self.orientation, "extremum")
        lon, lat = ext.shape(ishape).points[0]
        out = {}
        for layer in ("max_speed_contour", "outermost_contour"):
            r = self.shpc.get_reader(i_slice, self.orientation, layer)
            pts = np.asarray(r.shape(ishape).points)
            rad = r.record(ishape)["r_eq_area"]
            out[layer] = pts if (len(pts) > 2 and rad > 0) else None
            out[layer + "_r"] = rad if rad > 0 else np.nan
        ms = out["max_speed_contour"]
        if ms is None:
            ms = out["outermost_contour"]
        return dict(node=node, time=self.time(d), date_index=d, lon=lon, lat=lat,
                    contour=ms, outer=out["outermost_contour"],
                    r_maxspeed_km=out["max_speed_contour_r"],
                    r_outer_km=out["outermost_contour_r"],
                    speed=ext.record(ishape)["speed"], ssh=ext.record(ishape)["ssh"],
                    run=self.name)

    def nodes_at(self, date_index):
        i_slice = self.shpc.get_slice(date_index)
        if i_slice is None:
            return []
        self.shpc.get_ishape_last(i_slice, self.orientation)
        k = date_index - self.shpc.d_init[i_slice]
        e_max = int(self.shpc._e_max[i_slice][self.orientation][k])
        return [util_eddies.date_eddy_to_node(date_index, i, self.e_over)
                for i in range(1, e_max + 1)]

    def find_node(self, t, lon, lat, max_km=100.0, verbose=True):
        """Instantaneous eddy at time t whose contour contains (lon, lat), else the
        nearest center within max_km. Only eddies that are in the graph are kept."""
        d = self.date_index(t)
        cands = [self.eddy(n) for n in self.nodes_at(d) if n in self.node2seg]
        if not cands:
            raise ValueError(f"{self.name}: no {self.orientation} at {self.time(d)}")
        dist = np.array([haversine_km(lon, lat, e["lon"], e["lat"]) for e in cands])
        inside = [k for k, e in enumerate(cands)
                  for xy in (e["contour"], e["outer"])
                  if xy is not None and _points_in_poly(xy, lon, lat)]
        k = min(inside, key=lambda k: dist[k]) if inside else int(dist.argmin())
        if not inside and dist[k] > max_km:
            raise ValueError(f"{self.name}: no {self.orientation} within {max_km} km of "
                             f"({lon:.2f}, {lat:.2f}) at {self.time(d)} "
                             f"(nearest {dist[k]:.0f} km)")
        if verbose:
            how = "contour contains seed" if inside else "nearest center"
            print(f"  {self.name}: seed -> eddy at ({cands[k]['lon']:.2f}, "
                  f"{cands[k]['lat']:.2f}) on {cands[k]['time']:%Y-%m-%d %Hh}, "
                  f"{dist[k]:.0f} km away ({how})")
        return cands[k]["node"]


# =============================================================================
# Lineage
# =============================================================================
@dataclass
class Lineage:
    main: pd.DataFrame                      # one row per instantaneous eddy of the tracked eddy
    events: pd.DataFrame                    # merging / splitting events on the main track
    splits: list = field(default_factory=list)   # [DataFrame] eddies that split off
    merges: list = field(default_factory=list)   # [DataFrame] eddies that merged in
    seed: dict = field(default_factory=dict)


def _rows(run, segs, t_min=None, t_max=None, role="main", branch_id=None):
    rows = []
    for s in segs:
        for n in run.nodes(s):
            e = dict(run.eddy(n))
            if (t_min is not None and e["time"] < t_min) or \
               (t_max is not None and e["time"] > t_max):
                continue
            e.update(segment=s, role=role, branch=branch_id)
            rows.append(e)
    return pd.DataFrame(rows)


def _chain(run, seg, forward, t_lim, max_segs=200):
    """Follow the lowest-cost successor (forward) / predecessor (backward)."""
    segs = [seg]
    for _ in range(max_segs):
        nxt = run.succs(segs[-1]) if forward else run.preds(segs[-1])
        if not nxt:
            break
        cand = nxt[0][0]
        t_edge = run.time(util_eddies.node_to_date_eddy(
            run.nodes(cand)[0 if forward else -1], run.e_over, only_date=True))
        if t_lim is not None and ((forward and t_edge > t_lim) or (not forward and t_edge < t_lim)):
            break
        segs.append(cand)
    return segs if forward else segs[::-1]


def _event(run, kind, node, partners):
    e = run.eddy(node)
    return dict(time=e["time"], kind=kind, lon=e["lon"], lat=e["lat"], run=run.name,
                n_branches=len(partners) + 1, contour=e["contour"])


def trace_run(run, seed_node, backward=True, forward=True, t_min=None, t_max=None,
              branch_days=None):
    """Lineage of the eddy containing seed_node in one run.

    t_min / t_max  : clip the main track (None = as far as the graph goes).
    branch_days    : follow split / merged branches for at most this many days
                     (None = as long as they live, clipped to [t_min, t_max]).
    """
    seg0 = run.node2seg[seed_node]
    main, events, splits, merges = [seg0], [], [], []
    bl = pd.Timedelta(days=branch_days) if branch_days else None

    def add_split(parent, kept_child):
        others = [w for w, _ in run.succs(parent) if w != kept_child]
        if not others:
            return
        last = run.nodes(parent)[-1]
        events.append(_event(run, "splitting", last, others))
        t_ev = run.time(util_eddies.node_to_date_eddy(last, run.e_over, only_date=True))
        lim = min(filter(None, [t_max, t_ev + bl if bl else None]), default=None)
        for w in others:
            splits.append((w, _chain(run, w, True, lim), t_ev))

    def add_merge(child, kept_parent):
        others = [u for u, _ in run.preds(child) if u != kept_parent]
        if not others:
            return
        first = run.nodes(child)[0]
        events.append(_event(run, "merging", first, others))
        t_ev = run.time(util_eddies.node_to_date_eddy(first, run.e_over, only_date=True))
        lim = max(filter(None, [t_min, t_ev - bl if bl else None]), default=None)
        for u in others:
            merges.append((u, _chain(run, u, False, lim), t_ev))

    def seg_time(s, which):
        n = run.nodes(s)[0 if which == "first" else -1]
        return run.time(util_eddies.node_to_date_eddy(n, run.e_over, only_date=True))

    if backward:
        s = seg0
        while True:
            p = run.preds(s)
            if not p or (t_min is not None and seg_time(s, "first") <= t_min):
                break
            parent = p[0][0]
            add_merge(s, parent)
            add_split(parent, s)
            main.insert(0, parent)
            s = parent
    if forward:
        s = seg0
        while True:
            c = run.succs(s)
            if not c or (t_max is not None and seg_time(s, "last") >= t_max):
                break
            child = c[0][0]
            add_split(s, child)
            add_merge(child, s)
            main.append(child)
            s = child

    df = _rows(run, main, t_min, t_max)
    ev = pd.DataFrame(events, columns=["time", "kind", "lon", "lat", "run",
                                       "n_branches", "contour"])
    if len(ev):
        ev = ev.drop_duplicates(["time", "kind", "lon", "lat"]).sort_values("time")
        if t_min is not None:
            ev = ev[ev.time >= t_min]
        if t_max is not None:
            ev = ev[ev.time <= t_max]
    sp = [_rows(run, segs, t_ev, None, "split", f"{run.name}:{w}")
          for w, segs, t_ev in splits]
    mg = [_rows(run, segs, None, t_ev, "merge", f"{run.name}:{u}")
          for u, segs, t_ev in merges]
    return Lineage(df.sort_values("time").reset_index(drop=True), ev.reset_index(drop=True),
                   [b for b in sp if len(b)], [b for b in mg if len(b)])


def _clip(lin, t0, t1):
    """Keep only the part of a lineage in [t0, t1) (None = open)."""
    def ok(t):
        t = pd.Series(t) if isinstance(t, pd.Series) else t
        m = (t >= t0) if t0 is not None else (t == t)
        return m & ((t < t1) if t1 is not None else True)
    br = lambda bs: [b for b in bs if bool(ok(b.time.min())) or bool(ok(b.time.max()))]
    return Lineage(lin.main[ok(lin.main.time)], lin.events[ok(lin.events.time)]
                   if len(lin.events) else lin.events, br(lin.splits), br(lin.merges))


def load_runs(dataset, orientation):
    cfg = DATASETS[dataset]
    runs = []
    for d in cfg["dirs"]:
        try:
            runs.append(Run(DIR_TOEDDIES / d, orientation, cfg["epoch"], cfg["step_h"]))
        except FileNotFoundError as e:
            print(f"  skipping {d}: {e}")
    if not runs:
        raise RuntimeError(f"no usable TOEddies run for {dataset}")
    return sorted(runs, key=lambda r: r.t_first)


def track(dataset, orientation, seed_time, seed_lon, seed_lat, t_min=None, t_max=None,
          branch_days=None, stitch=True, verbose=True):
    """Lineage of the eddy at (seed_time, seed_lon, seed_lat), stitched across runs."""
    seed_time = pd.Timestamp(seed_time)
    t_min = pd.Timestamp(t_min) if t_min is not None else None
    t_max = pd.Timestamp(t_max) if t_max is not None else None
    runs = load_runs(dataset, orientation)
    if verbose:
        print("Runs:", *runs, sep="\n  ")
    covering = [k for k, r in enumerate(runs) if r.t_first <= seed_time <= r.t_last]
    if not covering:
        raise ValueError(f"seed time {seed_time} not covered by any run")
    # run in which the seed is the most central
    k0 = max(covering, key=lambda k: min(seed_time - runs[k].t_first,
                                         runs[k].t_last - seed_time))

    def go(k, t, lon, lat, backward, forward, lo, hi):
        node = runs[k].find_node(t, lon, lat, verbose=verbose)
        return trace_run(runs[k], node, backward, forward, lo, hi, branch_days)

    lin = go(k0, seed_time, seed_lon, seed_lat, True, True, t_min, t_max)
    parts = [(lin, None, None)]
    if stitch:
        # backwards into earlier runs
        k, cur = k0, lin
        while k > 0:
            prev = runs[k - 1]
            if prev.t_last <= runs[k].t_first:
                break
            t_sw = runs[k].t_first + (prev.t_last - runs[k].t_first) / 2
            t_sw = runs[k].time(runs[k].date_index(t_sw))
            if cur.main.empty or cur.main.time.min() > t_sw or (t_min and t_sw <= t_min):
                break
            row = cur.main.iloc[(cur.main.time - t_sw).abs().argmin()]
            if verbose:
                print(f"  stitching {runs[k].name} -> {prev.name} at {t_sw:%Y-%m-%d %Hh}")
            parts[-1] = (parts[-1][0], t_sw, parts[-1][2])
            try:
                cur = go(k - 1, row.time, row.lon, row.lat, True, False, t_min, row.time)
            except ValueError as e:
                print("  ", e)
                break
            parts.append((cur, None, t_sw))
            k -= 1
        # forwards into later runs
        k, cur = k0, lin
        while k < len(runs) - 1:
            nxt = runs[k + 1]
            if nxt.t_first >= runs[k].t_last:
                break
            t_sw = nxt.t_first + (runs[k].t_last - nxt.t_first) / 2
            t_sw = nxt.time(nxt.date_index(t_sw))
            if cur.main.empty or cur.main.time.max() < t_sw or (t_max and t_sw >= t_max):
                break
            row = cur.main.iloc[(cur.main.time - t_sw).abs().argmin()]
            if verbose:
                print(f"  stitching {runs[k].name} -> {nxt.name} at {t_sw:%Y-%m-%d %Hh}")
            i_cur = 0 if cur is lin else len(parts) - 1
            parts[i_cur] = (parts[i_cur][0], parts[i_cur][1], t_sw)
            try:
                cur = go(k + 1, row.time, row.lon, row.lat, False, True, row.time, t_max)
            except ValueError as e:
                print("  ", e)
                break
            parts.append((cur, t_sw, None))
            k += 1

    clipped = [_clip(l, a, b) for l, a, b in parts]
    main = pd.concat([c.main for c in clipped]).sort_values("time").reset_index(drop=True)
    evs = [c.events for c in clipped if len(c.events)]
    events = (pd.concat(evs).sort_values("time").reset_index(drop=True) if evs
              else clipped[0].events)
    return Lineage(main, events, sum((c.splits for c in clipped), []),
                   sum((c.merges for c in clipped), []),
                   seed=dict(time=seed_time, lon=seed_lon, lat=seed_lat,
                             dataset=dataset, orientation=orientation,
                             data_start=runs[0].t_first, data_end=runs[-1].t_last))
