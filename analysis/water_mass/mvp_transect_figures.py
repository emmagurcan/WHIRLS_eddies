"""
MVP transect figures (combines anom_ctd, anom_woa, epv_calc and anom_water_masses).

For each MVP transect, saves 9 figures in figures/mvp/Transect N/:
  1  Theta, S_A, sigma0 and SADCP cross-track velocity, with map
  2  relative vorticity zeta/f (SADCP)
  3  anomalies relative to WOA (July): Theta and S_A on isopycnals, sigma0 at constant pressure
  4  anomalies on isopycnals relative to the CTD background: Theta, S_A, isopycnal displacement
  5  Ertel PV with the WOA background (Delta EPV_z on isopycnals)
  6  Ertel PV with the CTD background (Delta EPV_z on isopycnals)
  7  water masses (Barabinot et al. 2025 definitions + Indian Ocean mode waters)
  8  theta-S diagram with the literature definitions of water_mass_literature.py
  9  sections: where each literature water mass is found (number of studies that agree)

Section orientation (the same for every transect and every calculation):
  +x points east (north for the mostly meridional transects 2-5), so east is always on the
  right of the plots. Tick labels are the distance from the first cast, so a westward
  transect reads max -> min from left to right.
  v_perp > 0 to the LEFT of +x (northward on east-west transects), so (x, v_perp, z) is
  right-handed and zeta = dv_perp/dx and EPV have their physical sign on every transect.

Usage (from analysis/):
  python mvp_transect_figures.py            all transects
  python mvp_transect_figures.py 12 13      only transects 12 and 13
"""
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from scipy import ndimage
from scipy.ndimage import uniform_filter1d
import cartopy.feature as cfeature

from tools_eg import *
from water_mass_literature import WM_LIT, WM_NAMES, WM_LIT_COLORS, lit_mask, ranged

OUT_BASE = BASE / "figures/mvp"
DIR_WOA = BASE / "data/woa"
DIR_CTD = BASE / "data/CTD"

ADT_LEVELS = np.arange(-0.3, 1.5 + 0.001, 0.02)       # map background
ADT_TICKS = [-0.3, 0, 0.3, 0.6, 0.9, 1.2, 1.5]
VEL_LIM = 0.8              # m/s, colour limit of v_perp
ZMAX = 1000                # m, bottom of the zeta and EPV figures

# zeta (epv_calc)
SMOOTH_X, SMOOTH_Z = 3, 3  # running means (bins) on v_perp before differentiating
PMV_FACTOR = 1             # 1 = Rudnick (Barabinot 2025), 2 = poor man's vorticity
RO_LIM = 0.3               # colour limit of zeta/f

# EPV (epv_calc)
RHO0 = 1025.0
DZ_EPV = 8.0               # m, depth bins of the EPV grid
SMOOTH_Z_B = 3             # running mean on b in depth (bins)
MAX_XGAP_KM = 10           # fill MVP gaps along x up to this size
REF_ZMIN = 10.0            # m, background levels shallower than this are ignored
MLD_DSIG = 0.03            # kg m^-3, background mixed-layer criterion
WOA_N2_SMOOTH_M = 150      # m, running mean on N^2_ref for WOA
EPV_X_MIN = 1e-13          # s^-3, R masked where |EPV_x| is below this
BETA = 40                  # eddy core: |Delta EPV_z| > BETA |EPV_x|
CORE_SMOOTH = 3
DEPV_MIN = 1e-11
CORE_ZMIN = 60             # m, automatic core seed searched below this depth
EPV_SCALE = 1e-10

# backgrounds
WOA_RADIUS_KM = 100        # WOA profiles within this distance of the eddy centre
REF_DP, REF_MIN_N = 1.0, 3 # CTD background: 1 dbar bins, levels with >= 3 casts


# ---------------------------------------------------------------------------
# Section geometry
# ---------------------------------------------------------------------------
def orient_section(sec):
    """+x towards the east (towards the north if the transect is mostly meridional).
    sec["x"] stays the distance from the first cast, with that sign."""
    ex = sec["xdir"]
    zonal = abs(ex[0]) >= abs(ex[1])
    sign = np.sign(ex[0] if zonal else ex[1]) or 1.0
    sec["x"], sec["xdir"] = sign * sec["x"], sign * ex
    sec["ends"] = ("W", "E") if zonal else ("S", "N")
    sec["vpos"] = "northward" if zonal else "westward"
    return sec


def finish_x_axis(sec, axs, map_ax=None):
    """Ticks = distance from the first cast; the label says which end is which."""
    w, e = sec["ends"]
    for ax in axs:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):g}"))
        if ax.get_xlabel():
            ax.set_xlabel(f"Distance from start of transect (km)      {w} ←   → {e}")
    if map_ax is not None:
        for t in map_ax.texts:
            if t.get_text().startswith("-"):
                t.set_text(t.get_text()[1:])


def v_perp_title(sec):
    return fr"(d) SADCP: v$_\perp$ {DS_KM:.0f} km bins ($>0$: {sec['vpos']})"


def mark_profiles(axs, x_km):
    """Small red dots just above the top of each panel, one per measured profile."""
    for ax in axs:
        ax.plot(x_km, np.full(len(x_km), 1.0), "o", ms=2.5, color="red", mec="none",
                transform=ax.get_xaxis_transform(), clip_on=False, zorder=10)


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


# ---------------------------------------------------------------------------
# Loading one transect
# ---------------------------------------------------------------------------
def merged_sadcp(sadcp_tr_all):
    """OS150 above Z_SPLIT + OS38 below, on the same 10-min time axis."""
    a = sadcp_window("OS150", sadcp=cascade_to_window_format(sadcp_tr_all["OS150"], pgood_min=50))
    b = sadcp_window("OS38", sadcp=cascade_to_window_format(sadcp_tr_all["OS38"], pgood_min=50))
    t = np.union1d(a["time"], b["time"])

    def on_t(d, key):
        out = np.full((len(t),) + d[key].shape[1:], np.nan)
        out[np.searchsorted(t, d["time"])] = d[key]
        return out

    lon_a, lat_a = on_t(a, "lon"), on_t(a, "lat")
    return dict(time=t,
                lon=np.where(np.isnan(lon_a), on_t(b, "lon"), lon_a),
                lat=np.where(np.isnan(lat_a), on_t(b, "lat"), lat_a),
                z=np.concatenate([a["z"], b["z"]]),
                u=np.concatenate([on_t(a, "u"), on_t(b, "u")], axis=1),
                v=np.concatenate([on_t(a, "v"), on_t(b, "v")], axis=1))


def load_transect(name, sadcp_full):
    """MVP section (oriented), merged SADCP, eddies crossed and map arguments."""
    files = sorted((DIR_MVP / name.lower()).glob("*.csv"))
    mvp = pd.concat([pd.read_csv(f).assign(CAST=k) for k, f in enumerate(files)], ignore_index=True)

    sadcp_tr_all = {inst: sadcp_for_transect(mvp, s) for inst, s in sadcp_full.items()}
    nav = sadcp_on_mvp_times(mvp, sadcp_tr_all["OS150"], max_gap="30min", keep=("LATITUDE", "LONGITUDE"))
    mvp["LON (°E)"] = mvp["LON (°E)"].fillna(pd.Series(nav["LONGITUDE"].values, index=mvp.index))
    mvp["LAT (°N)"] = mvp["LAT (°N)"].fillna(pd.Series(nav["LATITUDE"].values, index=mvp.index))
    adcp = merged_sadcp(sadcp_tr_all)

    date = pd.to_datetime(mvp["TIME"]).median().strftime("%Y-%m-%d")
    date_compact = date.replace("-", "")
    ssh = load_satellite(DIR_SSH / [f for f in daily_ssh_files if date in f or date_compact in f][0])
    eddy_files, eddy_source = get_eddy_files(date_compact, date)
    edds = eddies_on_transect(eddy_files, mvp["LON (°E)"].values, mvp["LAT (°N)"].values)

    casts = fill_cast_positions(split_by_cast_id(mvp, "CAST"))
    sec = build_section(casts, dp=1.0)
    sec["x"], sec["xdir"] = along_track_x(sec["lon"], sec["lat"])
    sec = orient_section(sec)

    # eddy used for the WOA background and the EPV core: the anticyclone closest to the track
    acs = [e for e in edds if not e["cyc"]]
    eddy = min(acs or edds, key=lambda e: e["dist_km"]) if edds else None

    crossed = ", ".join(("C" if e["cyc"] else "AC") + f" ({e['dist_km']:.0f} km)" for e in edds)
    print(f"{name}: {len(casts)} casts, {date}, {sec['ends'][0]}->{sec['ends'][1]} axis, "
          f"eddies crossed: {crossed or 'none'}")
    return dict(name=name, sec=sec, adcp=adcp, edds=edds, eddy=eddy, date=date,
                map_args=[ssh, eddy_files, name, eddy_source, ADT_TICKS, date, ADT_LEVELS])


# ---------------------------------------------------------------------------
# Backgrounds: WOA around the eddy, CTD mean of all casts
# ---------------------------------------------------------------------------
def load_woa():
    ds_t = nc.Dataset(sorted((DIR_WOA / "temp").glob("*.nc"))[0])
    ds_s = nc.Dataset(sorted((DIR_WOA / "salt").glob("*.nc"))[0])
    lon2d, lat2d = np.meshgrid(np.array(ds_t["lon"][:]), np.array(ds_t["lat"][:]))
    return dict(lon2d=lon2d, lat2d=lat2d, depth=np.array(ds_t["depth"][:]),
                temp=np.ma.filled(ds_t["t_an"][:], np.nan), sal=np.ma.filled(ds_s["s_an"][:], np.nan))


def woa_reference(woa, lon0, lat0):
    """Mean WOA profile within WOA_RADIUS_KM of (lon0, lat0) -> dict(p, CT, SA, sig0, lon, lat)."""
    dist = np.hypot((woa["lat2d"] - lat0) * 111.0,
                    (woa["lon2d"] - lon0) * 111.0 * np.cos(np.deg2rad(lat0)))
    mask = (dist <= WOA_RADIUS_KM) & np.isfinite(woa["temp"][0, 0]) & np.isfinite(woa["sal"][0, 0])
    t = np.squeeze(np.nanmean(woa["temp"][..., mask], axis=-1))
    s = np.squeeze(np.nanmean(woa["sal"][..., mask], axis=-1))
    lon_r, lat_r = woa["lon2d"][mask].mean(), woa["lat2d"][mask].mean()
    p = gsw.p_from_z(-woa["depth"], lat_r)
    SA = gsw.SA_from_SP(s, p, lon_r, lat_r)
    CT = gsw.CT_from_t(SA, t, p)
    ok = np.isfinite(SA) & np.isfinite(CT)
    return dict(p=p[ok], CT=CT[ok], SA=SA[ok], sig0=gsw.sigma0(SA, CT)[ok],
                lon=woa["lon2d"][mask], lat=woa["lat2d"][mask])


def read_cnv(path):
    names, rows, bad, in_data = [], [], -9.990e-29, False
    with open(path, errors="ignore") as f:
        for line in f:
            if in_data:
                parts = line.split()
                if len(parts) == len(names):
                    rows.append(parts)
            elif line.startswith("# name"):
                names.append(line.split("=", 1)[1].split(":")[0].strip())
            elif line.startswith("# bad_flag"):
                bad = float(line.split("=")[1])
            elif line.startswith("*END*"):
                in_data = True
    df = pd.DataFrame(np.asarray(rows, float), columns=names)
    return df.mask(np.isclose(df, bad))


def read_ctd_profile(path):
    """p (dbar), in-situ T, practical S of one cast, downcast only."""
    names = {"p": ["prDM", "prdM", "prSM", "pr", "PRES", "pres", "pressure", "PRESSURE"],
             "t": ["t090C", "t068C", "t190C", "TEMP", "temp", "temperature", "TEMPERATURE"],
             "s": ["sal00", "sal11", "PSAL", "psal", "salinity", "SALINITY"]}
    std = {"p": "sea_water_pressure", "t": "sea_water_temperature", "s": "sea_water_practical_salinity"}

    def pick(cols, key, fallback=None):
        low = {c.lower(): c for c in cols}
        return next((low[n.lower()] for n in names[key] if n.lower() in low), fallback)

    if path.suffix.lower() == ".nc":
        ds = xr.open_dataset(path)
        cols = list(ds.variables)
        by_std = {k: next((v for v in cols if ds[v].attrs.get("standard_name") == s), None) for k, s in std.items()}
        p, t, s = (np.asarray(ds[pick(cols, k, by_std[k])].values, float).squeeze() for k in "pts")
        ds.close()
    else:
        df = read_cnv(path) if path.suffix.lower() == ".cnv" else pd.read_csv(path)
        p, t, s = (df[pick(list(df.columns), k)].to_numpy(float) for k in "pts")
    imax = int(np.nanargmax(p))
    return p[:imax + 1], t[:imax + 1], s[:imax + 1]


def find_profile_file(name, md_nc):
    ext = (".nc", ".cnv", ".csv")
    stem = Path(name).stem
    hits = [f for f in DIR_CTD.rglob(f"{stem}*") if f.suffix.lower() in ext and f.resolve() != md_nc.resolve()]
    hits.sort(key=lambda f: (ext.index(f.suffix.lower()), f.stem != stem, len(f.name)))
    return hits[0] if hits else None


def ctd_reference():
    """Mean Theta / S_A of all CTD casts on 1 dbar bins -> dict(p, CT, SA, sig0, lon, lat)."""
    md_nc = DIR_CTD / "ctd_metadata.nc"
    casts = pd.read_csv(DIR_CTD / "ctd_toeddies_classification.csv")
    md = xr.open_dataset(md_nc)
    md_file = [v.decode().strip() if isinstance(v, bytes) else str(v).strip()
               for v in np.asarray(md["file"].values).ravel()]
    p_edges = np.arange(0, np.nanmax(md["max_pressure"].values) + REF_DP, REF_DP)
    md.close()
    p_mid = 0.5 * (p_edges[1:] + p_edges[:-1])

    def on_grid(p, val):
        idx = np.digitize(p, p_edges) - 1
        ok = (idx >= 0) & (idx < len(p_mid)) & np.isfinite(val)
        n = np.bincount(idx[ok], minlength=len(p_mid))
        with np.errstate(invalid="ignore"):
            return np.where(n > 0, np.bincount(idx[ok], val[ok], minlength=len(p_mid)) / n, np.nan)

    CT_all, SA_all, used = [], [], []
    for _, c in casts.iterrows():
        path = find_profile_file(md_file[c["cast"]], md_nc)
        if path is None:
            continue
        p, t, sp = read_ctd_profile(path)
        SA = gsw.SA_from_SP(sp, p, c["lon"], c["lat"])
        CT_all.append(on_grid(p, gsw.CT_from_t(SA, t, p)))
        SA_all.append(on_grid(p, SA))
        used.append(c)
    CT_all, SA_all, used = np.array(CT_all), np.array(SA_all), pd.DataFrame(used)
    keep = np.isfinite(CT_all).sum(axis=0) >= REF_MIN_N
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        CT, SA = np.nanmean(CT_all, axis=0), np.nanmean(SA_all, axis=0)
    print(f"CTD background: {len(used)} casts")
    return dict(p=p_mid[keep], CT=CT[keep], SA=SA[keep], sig0=gsw.sigma0(SA, CT)[keep],
                lon=used["lon"].values, lat=used["lat"].values, n_casts=len(used))


# ---------------------------------------------------------------------------
# Anomalies (anom_woa, anom_ctd)
# ---------------------------------------------------------------------------
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


def anom(ref, sec, var_name, pas=0.01):
    """Section anomaly relative to a background profile, same shape as sec["CT"].
    temp, sal, disp: on isopycnals (Barabinot et al. 2025, eq. 1-2), put back at the pressure
    of each isopycnal in each cast (disp = isopycnal depth - background depth, > 0 deeper).
    dens: sigma0 - sigma0_ref at constant pressure."""
    p, sig = sec["p"], sec["sig0"]
    if var_name == "dens":
        return sig - np.interp(p, ref["p"], ref["sig0"], left=np.nan, right=np.nan)[:, None]
    ref_profile, var_grid = {"temp": (ref["CT"], sec["CT"]),
                             "sal": (ref["SA"], sec["SA"]),
                             "disp": (ref["p"], np.broadcast_to(p[:, None], sig.shape))}[var_name]
    ind_rho = np.arange(22, 28, pas)
    ref_rho = interp_on_rho(ref["sig0"], ref_profile, ind_rho)
    out = np.full(sig.shape, np.nan)
    for i in range(sig.shape[1]):
        a_rho = interp_on_rho(sig[:, i], var_grid[:, i], ind_rho) - ref_rho
        p_rho = interp_on_rho(sig[:, i], p, ind_rho)
        good = np.isfinite(a_rho) & np.isfinite(p_rho)
        if good.sum() < 2:
            continue
        order = np.argsort(p_rho[good])
        d, idx = np.unique(p_rho[good][order], return_index=True)
        if d.size >= 2:
            out[:, i] = np.interp(p, d, a_rho[good][order][idx], left=np.nan, right=np.nan)
    return out


def plot_velocity_panel(fig, axc, sec, adcp):
    """SADCP cross-track velocity, 3 km bins (as in plot_section_with_map)."""
    s = along_track_x(np.asarray(adcp["lon"]), np.asarray(adcp["lat"]),
                      sec["lon"][0], sec["lat"][0], sec["xdir"])[0] / 1e3
    ex = sec["xdir"]
    vn = -adcp["u"] * ex[1] + adcp["v"] * ex[0]
    edges = np.arange(np.floor(np.nanmin(s)), np.ceil(np.nanmax(s)) + DS_KM, DS_KM)
    xg, V, _ = bin_columns(s, vn, edges)
    V = fill_z_gaps(V, adcp["z"])
    fill_lev = np.linspace(-VEL_LIM, VEL_LIM, N_FILL + 1)
    cf = axc.contourf(xg, adcp["z"], V.T, levels=fill_lev, cmap=cmocean.cm.balance, extend="both")
    nstep = int(np.ceil(VEL_LIM / C_LINE_STEP)) + 1
    cs = axc.contour(xg, adcp["z"], V.T, levels=np.round(np.arange(-nstep, nstep + 1) * C_LINE_STEP, 3),
                     colors="k", linewidths=0.6)
    axc.clabel(cs, cs.levels, fmt="%.1f", fontsize=6, inline=True, inline_spacing=2)
    cb = fig.colorbar(cf, ax=axc, pad=0.01)
    cb.set_ticks(fill_lev[1::3])
    cb.ax.set_yticklabels([f"{0.0 if abs(t) < 1e-9 else t:.2f}" for t in fill_lev[1::3]])
    cb.set_label(r"v$_\perp$ (m s$^{-1}$)")
    axc.set_title(v_perp_title(sec), loc="left", fontsize=11, pad=14)


def plot_anomalies(r, anoms, specs, title, ref, savepath, sig_min_ref=None):
    """(a)-(c) anomalies, (d) SADCP v_perp, small map on the right (orange = background points)."""
    sec, adcp = r["sec"], r["adcp"]
    lon, lat = sec["lon"], sec["lat"]
    xkm, p = sec["x"] / 1e3, sec["p"]
    sig_levels = np.arange(np.floor(np.nanmin(sec["sig0"]) * 4) / 4, np.nanmax(sec["sig0"]), 0.25)

    fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    outer = GridSpec(1, 2, figure=fig, width_ratios=[4, 1])
    gs_sec = outer[0].subgridspec(2, 2)
    gs_map = outer[1].subgridspec(2, 1, height_ratios=[1, 0.6])
    ax0 = fig.add_subplot(gs_sec[0, 0])
    axs = [ax0] + [fig.add_subplot(gs_sec[i, j], sharex=ax0) for i, j in [(0, 1), (1, 0), (1, 1)]]

    for ax, (key, lab, ttl) in zip(axs[:3], specs):
        vals = np.abs(anoms[key][np.isfinite(anoms[key])])
        lim = np.nanpercentile(vals, 99) if vals.size else 1.0
        pc = ax.pcolormesh(xkm, p, anoms[key], cmap=cmocean.cm.balance, shading="nearest", vmin=-lim, vmax=lim)
        fig.colorbar(pc, ax=ax, label=lab, pad=0.01)
        ax.set_title(ttl, loc="left", fontsize=11, pad=14)
    plot_velocity_panel(fig, axs[3], sec, adcp)
    if r["edds"]:
        mark_eddies_on_sections(axs, sec, adcp["z"], adcp["lon"], adcp["lat"], adcp["u"], adcp["v"], r["edds"])

    for ax in axs:
        cs = ax.contour(xkm, p, sec["sig0"], levels=sig_levels, colors="k", linewidths=0.5)
        ax.clabel(cs, cs.levels[::2], fontsize=7, fmt="%.2f")
        ax.plot(xkm, np.zeros_like(xkm), "|", color="k", ms=6, clip_on=False)
    if sig_min_ref is not None:            # lighter than this: no matching isopycnal in the background
        for ax in axs[:3]:
            ax.contour(xkm, p, sec["sig0"], levels=[sig_min_ref], colors="magenta", linewidths=1.2, linestyles="--")
    for ax in axs[:3]:
        ax.set_ylim(np.nanmax(p), 0)
        ax.set_ylabel("Pressure (dbar)")
    axs[3].set_ylim(np.nanmax(adcp["z"]), 0)
    axs[3].set_ylabel("Depth (m)")
    axs[0].set_xlim(np.nanmin(xkm), np.nanmax(xkm))
    for ax in axs[:2]:
        ax.tick_params(labelbottom=False)
    for ax in axs[2:]:
        ax.set_xlabel("x")

    # map
    crs = ccrs.PlateCarree()
    ext = [np.nanmin(lon) - 2.5, np.nanmax(lon) + 2.5, np.nanmin(lat) - 2.0, np.nanmax(lat) + 2.0]
    axm = fig.add_subplot(gs_map[0], projection=ccrs.Mercator(central_longitude=np.mean(ext[:2])))
    axm.set_extent(ext, crs=crs)
    draw_background(axm, crs, *r["map_args"])
    axm.add_feature(cfeature.LAND, facecolor="0.8", zorder=50)
    axm.coastlines(resolution="50m", linewidth=0.6, zorder=51)
    gl = axm.gridlines(crs=crs, draw_labels=True, linewidth=0.4, color="0.5", alpha=0.4, linestyle="--")
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 8}
    axm.plot(lon, lat, "-", color="k", lw=2, zorder=100, transform=crs,
             path_effects=[pe.Stroke(linewidth=4.5, foreground="white"), pe.Normal()])
    axm.scatter(ref["lon"], ref["lat"], s=12, color="orange", edgecolor="k", linewidth=0.3, zorder=101, transform=crs)
    for i, lab in [(0, "0 km"), (-1, f"{abs(xkm[-1]):.0f} km")]:
        axm.plot(lon[i], lat[i], marker="o" if i == 0 else "s", ms=6, mfc="yellow", mec="k", mew=1,
                 zorder=102, transform=crs)
        axm.annotate(lab, (lon[i], lat[i]), xytext=(6, 0), textcoords="offset points", fontsize=8,
                     fontweight="bold", va="center", zorder=103, xycoords=crs._as_mpl_transform(axm),
                     path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    axm.set_title(f"{sec['time'][0]:%Y-%m-%d %H:%M}\nto {sec['time'][-1]:%Y-%m-%d %H:%M} UTC", fontsize=9)

    finish_x_axis(sec, axs)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Vorticity and Ertel PV (epv_calc)
# ---------------------------------------------------------------------------
def vorticity(sec, adcp):
    """v_perp in 3 km bins along x, zeta = dv_perp/dx, f at each bin."""
    z = np.asarray(adcp["z"], float)
    s = along_track_x(np.asarray(adcp["lon"], float), np.asarray(adcp["lat"], float),
                      sec["lon"][0], sec["lat"][0], sec["xdir"])[0] / 1e3
    ex = sec["xdir"]
    v_perp = -adcp["u"] * ex[1] + adcp["v"] * ex[0]
    edges = np.arange(np.floor(np.nanmin(s)), np.ceil(np.nanmax(s)) + DS_KM, DS_KM)
    xg, V, _ = bin_columns(s, v_perp, edges)
    _, LATb, _ = bin_columns(s, np.asarray(adcp["lat"], float)[:, None], edges)
    V = nan_smooth(nan_smooth(fill_z_gaps(V, z), SMOOTH_X, axis=0), SMOOTH_Z, axis=1)
    zeta = PMV_FACTOR * np.gradient(V, xg * 1e3, axis=0)
    lat_b = np.where(np.isfinite(LATb[:, 0]), LATb[:, 0], np.nanmean(LATb[:, 0]))
    f = 2 * OMEGA * np.sin(np.deg2rad(lat_b))
    return dict(x=xg, z=z, v_perp=V, zeta=zeta, f=f, ro=zeta / f[:, None])


def ertel_pv(sec, vort):
    """EPV_z = (f + zeta) N^2 and EPV_x = -(dv_perp/dz)(db/dx) on (x bins, DZ_EPV depth bins)."""
    xg, z_a = vort["x"], vort["z"]
    edges_x = np.r_[xg - DS_KM / 2, xg[-1] + DS_KM / 2]
    depth_mvp = -gsw.z_from_p(sec["p"], np.nanmean(sec["lat"]))
    has_v = np.isfinite(vort["v_perp"]).any(axis=0)
    z_top = np.nanmin(z_a[has_v])
    z_bot = min(ZMAX, np.nanmax(z_a[has_v]), np.nanmax(depth_mvp[np.isfinite(sec["sig0"]).any(axis=1)]))
    z_edges = np.arange(z_top - DZ_EPV / 2, z_bot + DZ_EPV / 2, DZ_EPV)
    z = 0.5 * (z_edges[1:] + z_edges[:-1])

    V = np.full((len(xg), len(z)), np.nan)
    for i, row in enumerate(vort["v_perp"]):
        ok = np.isfinite(row)
        if ok.sum() >= 2:
            V[i] = np.interp(z, z_a[ok], row[ok], left=np.nan, right=np.nan)
    zeta = PMV_FACTOR * np.gradient(V, xg * 1e3, axis=0)
    dvdz = -np.gradient(V, z, axis=1)                                 # z up

    b_sec = -G * (sec["sig0"] + 1000.0 - RHO0) / RHO0
    _, Bx, _ = bin_columns(sec["x"] / 1e3, b_sec.T, edges_x)
    _, Bxz, _ = bin_columns(depth_mvp, Bx.T, z_edges)
    B = pd.DataFrame(Bxz.T).interpolate(axis=0, limit=int(MAX_XGAP_KM // DS_KM), limit_area="inside").values
    B = nan_smooth(nan_smooth(B, SMOOTH_X, axis=0), SMOOTH_Z_B, axis=1)
    N2 = -np.gradient(B, z, axis=1)
    dbdx = np.gradient(B, xg * 1e3, axis=0)

    EPV_z = (vort["f"][:, None] + zeta) * N2
    EPV_x = -dvdz * dbdx
    return dict(x=xg, z=z, f=vort["f"], EPV_z=EPV_z, EPV_x=EPV_x,
                sig0=-(B * RHO0 / G) + RHO0 - 1000.0, lat=np.nanmean(sec["lat"]))


def reference_N2(ref, lat0, dz=DZ_EPV):
    """Background sigma0 and N^2 on its own depth grid (CTD: binned; WOA: from its layers, smoothed)."""
    d = -gsw.z_from_p(ref["p"], lat0)
    b = -G * (ref["sig0"] + 1000.0 - RHO0) / RHO0
    o = np.argsort(d)
    d, b = d[o], b[o]
    d, b = d[d >= REF_ZMIN], b[d >= REF_ZMIN]
    edges = np.arange(REF_ZMIN, d.max() + dz, dz)
    zr = 0.5 * (edges[1:] + edges[:-1])
    if np.nanmedian(np.diff(d)) <= dz:                                 # dense (CTD)
        _, bb, _ = bin_columns(d, b[:, None], edges)
        bb = nan_smooth(bb[:, 0][None, :], SMOOTH_Z_B, axis=1)[0]
        N2 = nan_smooth(-np.gradient(bb, zr)[None, :], SMOOTH_Z_B, axis=1)[0]
    else:                                                              # coarse (WOA)
        bb = np.interp(zr, d, b, left=np.nan, right=np.nan)
        N2 = np.interp(zr, 0.5 * (d[1:] + d[:-1]), -np.diff(b) / np.diff(d), left=np.nan, right=np.nan)
        N2 = nan_smooth(N2[None, :], max(1, int(round(WOA_N2_SMOOTH_M / dz))), axis=1)[0]
    return dict(sig0=-(bb * RHO0 / G) + RHO0 - 1000.0, N2=N2)


def epv_anomaly(epv, ref):
    """Delta EPV_z = EPV_z - f N^2_ref, with N^2_ref taken ON THE SAME ISOPYCNAL (Barabinot eqs. 11, 14)."""
    prof = reference_N2(ref, epv["lat"])
    ok = np.isfinite(prof["sig0"]) & np.isfinite(prof["N2"])
    s, i = np.unique(np.maximum.accumulate(prof["sig0"][ok]), return_index=True)
    N2r = np.interp(epv["sig0"], s, prof["N2"][ok][i], left=np.nan, right=np.nan)
    dEPV_z = epv["EPV_z"] - epv["f"][:, None] * N2r
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(np.abs(epv["EPV_x"]) > EPV_X_MIN, dEPV_z / epv["EPV_x"], np.nan)
    return dict(dEPV_z=dEPV_z, R=R)


def eddy_core(epv, an, eddy, to_sn):
    """Connected region around the eddy centre where |Delta EPV_z| > BETA |EPV_x| (Barabinot sect. 4.5)."""
    if eddy is None:
        return None
    sgn = -1.0 if eddy["cyc"] else 1.0                     # SH: AC -> Delta EPV_z > 0
    dE = sgn * nan_smooth(nan_smooth(an["dEPV_z"], CORE_SMOOTH, 0), CORE_SMOOTH, 1)
    Ex = nan_smooth(nan_smooth(np.abs(epv["EPV_x"]), CORE_SMOOTH, 0), CORE_SMOOTH, 1)
    crit = np.isfinite(dE) & np.isfinite(Ex) & (dE > BETA * Ex) & (dE > DEPV_MIN)
    ix = int(np.argmin(np.abs(epv["x"] - float(to_sn(*eddy["center"])[0]))))
    col = np.where(epv["z"] >= CORE_ZMIN, sgn * an["dEPV_z"][ix], np.nan)
    if not np.isfinite(col).any():
        return None
    lab, _ = ndimage.label(crit, structure=np.ones((3, 3)))
    seed = lab[ix, int(np.nanargmax(col))]
    return ndimage.binary_fill_holes(lab == seed) if seed else None


def plot_zeta(r, vort, savepath):
    sec = r["sec"]
    fig, ax = plt.subplots(figsize=(11, 4.5), constrained_layout=True)
    pc = ax.contourf(vort["x"], vort["z"], vort["ro"].T, cmap="RdBu_r", extend="both",
                     levels=np.linspace(-RO_LIM, RO_LIM, int(RO_LIM * 200) + 1))
    sig_lev = np.arange(np.floor(np.nanmin(sec["sig0"]) * 4) / 4, np.nanmax(sec["sig0"]), 0.25)
    cs = ax.contour(sec["x"] / 1e3, sec["p"], sec["sig0"], levels=sig_lev, colors="k", linewidths=0.6, alpha=0.8)
    ax.clabel(cs, fmt="%.2f", fontsize=6)
    ax.set_ylim(ZMAX, 0)
    ax.set_xlim(np.nanmin(sec["x"]) / 1e3, np.nanmax(sec["x"]) / 1e3)
    ax.set_xlabel("x")
    ax.set_ylabel("Depth (m)")
    fig.colorbar(pc, ax=ax, pad=0.01).set_label(r"$\zeta / f$")
    ax.set_title(rf"MVP {r['name']} ({r['date']}): $\zeta/f$  ($<0$: anticyclonic)", loc="left", fontsize=10)
    finish_x_axis(sec, [ax])
    fig.savefig(savepath, dpi=300)
    plt.close(fig)


def plot_epv(r, epv, an, core, title, savepath):
    """(a) EPV_z  (b) EPV_x  (c) Delta EPV_z  (d) |Delta EPV_z / EPV_x|, eddy core in red."""
    x, z = epv["x"], epv["z"]
    fig, axs = plt.subplots(2, 2, figsize=(16, 8), sharex=True, sharey=True, constrained_layout=True)
    axs = axs.ravel()
    for ax, (F, ttl) in zip(axs[:3], [(epv["EPV_z"], r"(a) $EPV_z$"), (epv["EPV_x"], r"(b) $EPV_x$"),
                                      (an["dEPV_z"], r"(c) $\Delta EPV_z$ (on isopycnals)")]):
        Fs = F / EPV_SCALE
        lim = np.nanpercentile(np.abs(Fs), 98)
        pc = ax.contourf(x, z, Fs.T, levels=np.linspace(-lim, lim, 41), cmap=cmocean.cm.balance, extend="both")
        fig.colorbar(pc, ax=ax, pad=0.01, ticks=np.linspace(-lim, lim, 5), format="%.2g").set_label(r"$10^{-10}$ s$^{-3}$")
        ax.set_title(ttl, loc="left", fontsize=11)
    pc = axs[3].contourf(x, z, np.abs(an["R"]).T, levels=np.linspace(0, 2 * BETA, 41), cmap=cmocean.cm.amp, extend="max")
    cb = fig.colorbar(pc, ax=axs[3], pad=0.01, ticks=np.linspace(0, 2 * BETA, 5))
    cb.set_label(r"$|\Delta EPV_z / EPV_x|$")
    cb.ax.axhline(BETA, color="k", lw=1.5)
    axs[3].set_title(r"(d) $|\Delta EPV_z / EPV_x|$", loc="left", fontsize=11)

    sig_lev = np.arange(np.floor(np.nanmin(epv["sig0"]) * 4) / 4, np.nanmax(epv["sig0"]), 0.25)
    for ax in axs:
        cs = ax.contour(x, z, epv["sig0"].T, levels=sig_lev, colors="k", linewidths=0.5)
        ax.clabel(cs, cs.levels, fontsize=6, fmt="%.2f")
        if core is not None and core.any():
            ax.contour(x, z, core.T.astype(float), levels=[0.5], colors="red", linewidths=1.8,
                       path_effects=[pe.Stroke(linewidth=3.2, foreground="white"), pe.Normal()])
        ax.set_ylim(min(ZMAX, z[-1]), 0)
    for ax in axs[2:]:
        ax.set_xlabel("x")
    for ax in axs[::2]:
        ax.set_ylabel("Depth (m)")
    if core is not None and core.any():
        title += rf"   (red: eddy core, $|\Delta EPV_z / EPV_x| > {BETA}$)"
    finish_x_axis(r["sec"], axs)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.savefig(savepath, dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Water masses (anom_water_masses)
# ---------------------------------------------------------------------------
def _mode_water_box(pt, dpt, sp, dsp, sig, dsig, lon, lat):
    """Theta / S / sigma_theta box (mean +/- 1 std) converted to Theta / S_A."""
    SA = gsw.SA_from_SP(np.array([sp - dsp, sp + dsp]), 0, lon, lat)
    CT = gsw.CT_from_pt(SA.mean(), np.array([pt - dpt, pt + dpt]))
    return dict(CT=tuple(CT), SA=tuple(SA), sig=(sig - dsig, sig + dsig))


ESACW = dict(CT=(9.44, 13.60), SA=(34.90, 35.40), sig=(26.5, 26.93))
IOSTMW = _mode_water_box(16.54, 0.49, 35.31, 0.04, 26.00, 0.10, lon=37.5, lat=-32.5)
IOSTMW2 = _mode_water_box(12.84, 0.57, 35.17, 0.11, 26.57, 0.04, lon=70.0, lat=-39.5)
AAIW_SIG_MIN = 26.9
WM_SHORT = ["SMW", "SACW", "IOSTMW", "ESACW", "IOSTMW2", "AAIW"]
WM_LABEL = ["SMW\n(σ0 < 25.6)",
            "SACW\n(σ0 25.6–26.5)",
            "IOSTMW\n(θ 16.54±0.49 °C,\nS 35.31±0.04,\nσθ 26.0±0.1)",
            "ESACW\n(Θ 9.44–13.60 °C,\nS_A 34.90–35.40,\nσ0 26.5–26.93)",
            "IOSTMW2\n(θ 12.84±0.57 °C,\nS 35.17±0.11,\nσθ 26.57±0.04)",
            "AAIW\n(σ0 ≥ 26.9)"]
WM_COLORS = ["#fde68a", "#f59e0b", "#0d9488", "#ea580c", "#2563eb", "#7c3aed"]
WM = {name: k for k, name in enumerate(WM_SHORT)}


def _in_box(box, CT, SA, sig0):
    return ((sig0 >= box["sig"][0]) & (sig0 <= box["sig"][1]) & (CT >= box["CT"][0])
            & (CT <= box["CT"][1]) & (SA >= box["SA"][0]) & (SA <= box["SA"][1]))


def classify_water_masses(CT, SA, sig0):
    """Index in WM_SHORT, NaN where no definition fits. Later assignments win."""
    cls = np.full(np.shape(sig0), np.nan)
    in_esacw = _in_box(ESACW, CT, SA, sig0)
    cls[sig0 < 25.6] = WM["SMW"]
    cls[(sig0 >= 25.6) & (sig0 < 26.5)] = WM["SACW"]
    cls[(sig0 >= AAIW_SIG_MIN) & ~in_esacw] = WM["AAIW"]
    cls[in_esacw] = WM["ESACW"]
    cls[_in_box(IOSTMW, CT, SA, sig0)] = WM["IOSTMW"]
    cls[_in_box(IOSTMW2, CT, SA, sig0)] = WM["IOSTMW2"]
    cls[~(np.isfinite(CT) & np.isfinite(SA) & np.isfinite(sig0))] = np.nan
    return cls


def plot_water_masses(r, savepath):
    sec = r["sec"]
    xkm = sec["x"] / 1e3
    wm = classify_water_masses(sec["CT"], sec["SA"], sec["sig0"])
    fig, ax = plt.subplots(figsize=(13, 4.2), constrained_layout=True)
    pc = ax.pcolormesh(xkm, sec["p"], wm, cmap=ListedColormap(WM_COLORS),
                       norm=BoundaryNorm(np.arange(-0.5, len(WM_SHORT)), len(WM_SHORT)), shading="nearest")
    sig_lev = np.arange(np.floor(np.nanmin(sec["sig0"]) * 4) / 4, np.nanmax(sec["sig0"]), 0.25)
    cs = ax.contour(xkm, sec["p"], sec["sig0"], levels=sig_lev, colors="k", linewidths=0.5)
    ax.clabel(cs, cs.levels[::2], fontsize=7, fmt="%.2f")
    to_sn = make_to_sn(sec)
    for e in r["edds"]:                    # TOEddies outer contours crossing the section
        for xv in crossings_on_line(e["outer"], to_sn):
            ax.axvline(xv, color="tab:blue" if e["cyc"] else "tab:red", ls="--", lw=1.4)
    ax.set_ylim(np.nanmax(sec["p"]), 0)
    ax.set_xlim(np.nanmin(xkm), np.nanmax(xkm))
    ax.set_ylabel("Pressure (dbar)")
    ax.set_xlabel("x")
    ax.set_title(f"MVP {r['name']} ({r['date']}): water masses  (dashed: TOEddies outer contours, "
                 "red AC, blue C)", loc="left", fontsize=11)
    cb = fig.colorbar(pc, ax=ax, ticks=range(len(WM_SHORT)), pad=0.01)
    cb.ax.set_yticklabels(WM_LABEL, fontsize=8)
    finish_x_axis(sec, [ax])
    fig.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Literature water-mass definitions (water_mass_literature.py), in theta / S_P
# ---------------------------------------------------------------------------
def theta_sp(sec):
    """Potential temperature (°C) and practical salinity of the section, as in the literature tables."""
    P = np.broadcast_to(sec["p"][:, None], sec["SA"].shape)
    SP = gsw.SP_from_SA(sec["SA"], P, sec["lon"][None, :], sec["lat"][None, :])
    return gsw.pt_from_CT(sec["SA"], sec["CT"]), SP


def plot_lit_ts(r, savepath):
    """theta-S diagram of the section with every literature box (one colour per water mass,
    one line style per study) and the single-value definitions as stars."""
    sec = r["sec"]
    pt, SP = theta_sp(sec)
    P = np.broadcast_to(sec["p"][:, None], pt.shape)
    ok = np.isfinite(pt) & np.isfinite(SP)
    step = max(1, ok.sum() // 60000)
    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)
    sc = ax.scatter(SP[ok][::step], pt[ok][::step], c=P[ok][::step], s=1, cmap=cmocean.cm.deep, rasterized=True)
    fig.colorbar(sc, ax=ax, pad=0.01).set_label("Pressure (dbar)")
    # axis limits from the 0.5-99.5 percentiles, so that one bad cast does not stretch the plot
    s_lim = (np.nanpercentile(SP, 0.5) - 0.15, np.nanpercentile(SP, 99.5) + 0.15)
    t_lim = (np.nanpercentile(pt, 0.5) - 1.0, np.nanpercentile(pt, 99.5) + 1.5)

    # sigma0 contours in (S_P, theta) at the mean section position
    Sg, Tg = np.meshgrid(np.linspace(*s_lim, 200), np.linspace(*t_lim, 200))
    SAg = gsw.SA_from_SP(Sg, 0, np.nanmean(sec["lon"]), np.nanmean(sec["lat"]))
    cs = ax.contour(Sg, Tg, gsw.sigma0(SAg, gsw.CT_from_pt(SAg, Tg)), levels=np.arange(22, 28.01, 0.25),
                    colors="0.6", linewidths=0.5)
    ax.clabel(cs, cs.levels[::2], fontsize=7, fmt="%.2f")

    styles = ["-", "--", ":", "-."]
    handles, k_style = [], {}
    for e in WM_LIT:
        col = WM_LIT_COLORS[e["wm"]]
        lab = f"{e['wm']}: {e['study']}"
        if e["T"] is not None and e["S"] is not None:
            ls = styles[k_style.setdefault(e["wm"], 0) % len(styles)]
            k_style[e["wm"]] += 1
            s0, s1 = (v if v is not None else lim for v, lim in zip(e["S"], (-1e3, 1e3)))
            t0, t1 = (v if v is not None else lim for v, lim in zip(e["T"], (-1e3, 1e3)))
            h = mpatches.Rectangle((s0, t0), s1 - s0, t1 - t0, fill=False, ec=col, lw=1.8, ls=ls, label=lab)
            ax.add_patch(h)
            handles.append(h)
        elif e["T0"] is not None:
            handles.append(ax.plot(e["S0"], e["T0"], "*", ms=13, mfc=col, mec="k", mew=0.6,
                                   ls="none", label=lab + " (core)", clip_on=True)[0])
    ax.set_xlim(*s_lim)
    ax.set_ylim(*t_lim)
    ax.set_xlabel(r"Practical salinity $S_P$")
    ax.set_ylabel(r"Potential temperature $\theta$ (°C)")
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.12, 1), fontsize=8, frameon=False,
              title="Literature definitions\n(boxes / stars outside the axes are off the data range)",
              title_fontsize=8)
    ax.set_title(f"MVP {r['name']} ({r['date']}): $\\theta$-S with literature water-mass definitions",
                 loc="left", fontsize=11)
    fig.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_lit_sections(r, savepath):
    """One panel per literature water mass: number of studies whose (theta, S) box contains each point."""
    sec = r["sec"]
    pt, SP = theta_sp(sec)
    xkm, p = sec["x"] / 1e3, sec["p"]
    wms = [wm for wm in WM_NAMES if ranged(wm)]
    ncol = 2
    nrow = int(np.ceil(len(wms) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 2.8 * nrow), sharex=True, sharey=True,
                            constrained_layout=True, squeeze=False)
    axs = axs.ravel()
    sig_lev = np.arange(np.floor(np.nanmin(sec["sig0"]) * 4) / 4, np.nanmax(sec["sig0"]), 0.25)
    for ax, wm in zip(axs, wms):
        entries = ranged(wm)
        n = sum(lit_mask(e, pt, SP).astype(int) for e in entries).astype(float)
        n[~(np.isfinite(pt) & np.isfinite(SP))] = np.nan
        frac = np.mean(n[np.isfinite(n)] > 0) * 100
        ax.pcolormesh(xkm, p, np.where(np.isfinite(pt), 0.0, np.nan), cmap=ListedColormap(["0.92"]),
                      shading="nearest")
        cmap = ListedColormap([WM_LIT_COLORS[wm]]) if len(entries) == 1 else \
            ListedColormap(plt.matplotlib.colors.to_rgba_array([WM_LIT_COLORS[wm]] * len(entries))
                           * np.c_[np.ones((len(entries), 3)), np.linspace(0.35, 1, len(entries))])
        pc = ax.pcolormesh(xkm, p, np.where(n > 0, n, np.nan), cmap=cmap, shading="nearest",
                           norm=BoundaryNorm(np.arange(0.5, len(entries) + 1), len(entries)))
        cb = fig.colorbar(pc, ax=ax, pad=0.01, ticks=range(1, len(entries) + 1))
        cb.set_label("studies", fontsize=8)
        cs = ax.contour(xkm, p, sec["sig0"], levels=sig_lev, colors="k", linewidths=0.4)
        ax.clabel(cs, cs.levels[::4], fontsize=6, fmt="%.2f")
        studies = "; ".join(e["study"] for e in entries)
        ax.set_title(f"{wm} ({WM_NAMES[wm]}): {frac:.0f}% of the section\n{studies}",
                     loc="left", fontsize=9)
        ax.set_ylim(np.nanmax(p), 0)
        ax.set_xlim(np.nanmin(xkm), np.nanmax(xkm))
    for ax in axs[len(wms):]:
        ax.set_visible(False)
    for ax in axs[::ncol]:
        ax.set_ylabel("Pressure (dbar)")
    for ax in axs[len(wms) - ncol:len(wms)]:
        ax.set_xlabel("x")
        ax.tick_params(labelbottom=True)
    finish_x_axis(sec, axs[:len(wms)])
    fig.suptitle(f"MVP {r['name']} ({r['date']}): literature water-mass definitions ($\\theta$, $S_P$ boxes; "
                 "darker = more studies agree)", fontsize=12, fontweight="bold")
    fig.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# One transect -> 9 figures
# ---------------------------------------------------------------------------
def make_figures(r, ctd_ref, woa):
    name, sec, adcp = r["name"], r["sec"], r["adcp"]
    out = OUT_BASE / name
    out.mkdir(parents=True, exist_ok=True)

    # 1. hydrography + SADCP
    fig, axs = plot_section_with_map(sec, adcp, map_draw=draw_background, map_args=r["map_args"],
                                     title=f"MVP {name}", edd=r["edds"], flip_x=False, vel_lim=VEL_LIM)
    axs[4].set_title(v_perp_title(sec), loc="left", fontsize=11, pad=14)
    mark_profiles(axs[1:4], sec["x"] / 1e3)                       # MVP casts
    s_adcp = along_track_x(np.asarray(adcp["lon"]), np.asarray(adcp["lat"]),
                           sec["lon"][0], sec["lat"][0], sec["xdir"])[0] / 1e3
    mark_profiles([axs[4]], s_adcp[np.isfinite(adcp["v"]).any(axis=1)])   # 10-min SADCP ensembles
    finish_x_axis(sec, axs[1:], map_ax=axs[0])
    fig.savefig(out / f"{name}_1_sections_hydro_adcp.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 2. zeta / f
    vort = vorticity(sec, adcp)
    plot_zeta(r, vort, out / f"{name}_2_zeta_over_f.png")

    # 3. anomalies relative to WOA, around the eddy (or the transect centre if no eddy is crossed)
    lon0, lat0 = r["eddy"]["center"] if r["eddy"] else (np.nanmean(sec["lon"]), np.nanmean(sec["lat"]))
    woa_ref = woa_reference(woa, lon0, lat0)
    anoms = dict(CT=anom(woa_ref, sec, "temp"), SA=anom(woa_ref, sec, "sal"), sig0=anom(woa_ref, sec, "dens"))
    plot_anomalies(r, anoms,
                   [("CT", r"$\Delta\Theta$ (°C)", r"(a) $\Theta$ anomaly on isopycnals"),
                    ("SA", r"$\Delta S_A$ (g kg$^{-1}$)", r"(b) $S_A$ anomaly on isopycnals"),
                    ("sig0", r"$\Delta\sigma_0$ (kg m$^{-3}$)", r"(c) $\sigma_0$ anomaly at constant pressure")],
                   f"MVP {name}: anomalies relative to WOA (July)", woa_ref,
                   out / f"{name}_3_anomalies_woa.png", sig_min_ref=np.nanmin(woa_ref["sig0"]))

    # 4. anomalies on isopycnals relative to the CTD background
    anoms = dict(CT=anom(ctd_ref, sec, "temp"), SA=anom(ctd_ref, sec, "sal"), disp=anom(ctd_ref, sec, "disp"))
    plot_anomalies(r, anoms,
                   [("CT", r"$\Delta\Theta$ (°C)", r"(a) $\Theta$ anomaly on isopycnals"),
                    ("SA", r"$\Delta S_A$ (g kg$^{-1}$)", r"(b) $S_A$ anomaly on isopycnals"),
                    ("disp", r"$\Delta p$ (dbar)", r"(c) Isopycnal displacement ($>0$: deeper than background)")],
                   f"MVP {name}: anomalies on isopycnals relative to {ctd_ref['n_casts']} CTD casts", ctd_ref,
                   out / f"{name}_4_anomalies_ctd.png")

    # 5-6. Ertel PV, Delta EPV_z on isopycnals
    epv = ertel_pv(sec, vort)
    to_sn = make_to_sn(sec)
    for k, (ref, txt, tag) in enumerate([(woa_ref, "WOA July", "woa"),
                                         (ctd_ref, f"{ctd_ref['n_casts']} CTD casts", "ctd")]):
        an = epv_anomaly(epv, ref)
        plot_epv(r, epv, an, eddy_core(epv, an, r["eddy"], to_sn),
                 f"MVP {name}: EPV background = {txt}", out / f"{name}_{5 + k}_epv_{tag}.png")

    # 7. water masses
    plot_water_masses(r, out / f"{name}_7_water_masses.png")

    # 8-9. literature water-mass definitions (water_mass_literature.py)
    plot_lit_ts(r, out / f"{name}_8_lit_theta_S.png")
    plot_lit_sections(r, out / f"{name}_9_lit_water_mass_sections.png")
    print(f"  -> 9 figures in {out}")


if __name__ == "__main__":
    numbers = sys.argv[1:] or sorted((d.name.split()[-1] for d in DIR_MVP.glob("transect *")), key=int)
    transects = [f"Transect {n}" for n in numbers]

    sadcp_full = {}
    for inst in ["OS150", "OS38"]:
        s = cascade_add_time(xr.open_dataset(SADCP_FILES[inst])).sortby("time")
        sadcp_full[inst] = s.isel(time=~s.indexes["time"].duplicated())
    ctd_ref = ctd_reference()
    woa = load_woa()

    failed = {}
    for name in transects:
        try:
            make_figures(load_transect(name, sadcp_full), ctd_ref, woa)
        except Exception as err:
            failed[name] = f"{type(err).__name__}: {err}"
            print(f"  !! {name} failed: {failed[name]}")
            plt.close("all")
    print(f"\ndone: {len(transects) - len(failed)} of {len(transects)} transects")
    for name, msg in failed.items():
        print(f"  {name}: {msg}")
