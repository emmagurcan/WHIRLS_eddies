"""
CTD background (reference) profile for the WHIRLS anomaly / EPV calculations.

Computed ONCE from the CTD casts and saved to data/CTD/ctd_reference_profile.nc, so the
notebooks (epv_calc, anom_ctd, ...) only have to load it with load_ctd_reference().

Run:  python make_ctd_reference_profile.py            (all casts)
      python make_ctd_reference_profile.py --outside  (only casts outside every TOEddies contour,
                                                       saved as ctd_reference_profile_outside.nc)

The file contains
  dims  pressure (1 dbar bins, levels with >= REF_MIN_N casts), cast
  CT, SA, sig0, CT_std, SA_std, n   (pressure)   mean background and spread
  CT_cast, SA_cast                  (cast, pressure)   every cast on the same grid
  station, lon, lat, time, file     (cast)       casts actually used
"""
import argparse
import warnings
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import xarray as xr

BASE = Path("/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS")
DIR_CTD = BASE / "data/CTD"
MD_CTD_NC = DIR_CTD / "ctd_metadata.nc"
CLASS_CSV = DIR_CTD / "ctd_toeddies_classification.csv"
OUT_ALL = DIR_CTD / "ctd_reference_profile.nc"
OUT_OUTSIDE = DIR_CTD / "ctd_reference_profile_outside.nc"

PROFILE_EXT = (".nc", ".cnv", ".csv")   # preferred profile format first
REF_DP = 1.0                            # dbar, same bins as the MVP section
REF_MIN_N = 3                           # keep only pressure levels with >= this many casts

NAMES = {"p": ["prDM", "prdM", "prSM", "pr", "PRES", "pres", "pressure", "PRESSURE"],
         "t": ["t090C", "t068C", "t190C", "TEMP", "temp", "temperature", "TEMPERATURE"],
         "s": ["sal00", "sal11", "PSAL", "psal", "salinity", "SALINITY"]}
STD_NAMES = {"p": "sea_water_pressure", "t": "sea_water_temperature",
             "s": "sea_water_practical_salinity"}


def _txt(v):
    return v.decode().strip() if isinstance(v, bytes) else str(v).strip()


def find_profile_file(name):
    """Profile file of one cast under DIR_CTD, from the metadata 'file' entry."""
    stem = Path(name).stem
    excluded = {MD_CTD_NC.resolve(), OUT_ALL.resolve(), OUT_OUTSIDE.resolve()}
    hits = [f for f in DIR_CTD.rglob(f"{stem}*") if f.suffix.lower() in PROFILE_EXT
            and f.resolve() not in excluded]
    if not hits:
        return None
    hits.sort(key=lambda f: (PROFILE_EXT.index(f.suffix.lower()), f.stem != stem, len(f.name)))
    return hits[0]


def _pick(cols, key, std=None):
    for n in NAMES[key]:
        if n in cols:
            return n
    low = {c.lower(): c for c in cols}
    for n in NAMES[key]:
        if n.lower() in low:
            return low[n.lower()]
    return std


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
    """-> p (dbar), in-situ T, practical S of one cast (downcast only)."""
    suf = path.suffix.lower()
    if suf == ".nc":
        ds = xr.open_dataset(path)
        cols = list(ds.variables)
        by_std = {k: next((v for v in cols if ds[v].attrs.get("standard_name") == s), None)
                  for k, s in STD_NAMES.items()}
        names = {k: _pick(cols, k, by_std[k]) for k in "pts"}
        if None in names.values():
            raise KeyError(f"{path.name}: no pressure/temperature/salinity among {cols}")
        p, t, s = (np.asarray(ds[names[k]].values, float).squeeze() for k in "pts")
        ds.close()
    else:
        df = read_cnv(path) if suf == ".cnv" else pd.read_csv(path)
        names = {k: _pick(list(df.columns), k) for k in "pts"}
        if None in names.values():
            raise KeyError(f"{path.name}: no pressure/temperature/salinity among {list(df.columns)}")
        p, t, s = (df[names[k]].to_numpy(float) for k in "pts")
    if p.ndim != 1:
        raise ValueError(f"{path.name}: pressure is {p.ndim}-D, expected one profile per file")
    imax = int(np.nanargmax(p))                    # downcast: surface -> deepest point
    return p[:imax + 1], t[:imax + 1], s[:imax + 1], names


def bin_on_grid(p, val, p_edges):
    idx = np.digitize(p, p_edges) - 1
    ok = (idx >= 0) & (idx < len(p_edges) - 1) & np.isfinite(val)
    s = np.bincount(idx[ok], val[ok], minlength=len(p_edges) - 1)
    n = np.bincount(idx[ok], minlength=len(p_edges) - 1)
    with np.errstate(invalid="ignore"):
        return np.where(n > 0, s / n, np.nan)


def make_reference(use_all_casts=True):
    cls = pd.read_csv(CLASS_CSV, parse_dates=["time"])
    ref_casts = (cls if use_all_casts else cls[cls["outside_eddies"]]).reset_index(drop=True)
    md = xr.open_dataset(MD_CTD_NC)
    md_file = [_txt(v) for v in np.asarray(md["file"].values).ravel()]
    md_pmax = np.asarray(md["max_pressure"].values, float).ravel()
    md.close()
    print(f"{len(ref_casts)} casts used for the background "
          f"({'all casts' if use_all_casts else 'outside every TOEddies contour'})")

    p_edges = np.arange(0, np.nanmax(md_pmax) + REF_DP, REF_DP)
    p_ref = 0.5 * (p_edges[1:] + p_edges[:-1])
    CT_list, SA_list, used, files = [], [], [], []
    for _, c in ref_casts.iterrows():
        path = find_profile_file(md_file[c["cast"]])
        if path is None:
            print(f"  ! station {c['station']}: no profile file for '{md_file[c['cast']]}' under {DIR_CTD}")
            continue
        p, t, sp, names = read_ctd_profile(path)
        SA = gsw.SA_from_SP(sp, p, c["lon"], c["lat"])
        CT = gsw.CT_from_t(SA, t, p)
        CT_list.append(bin_on_grid(p, CT, p_edges))
        SA_list.append(bin_on_grid(p, SA, p_edges))
        used.append(c)
        files.append(path.name)
        print(f"  station {c['station']:>4}: {path.name} ({names['p']}, {names['t']}, {names['s']}), "
              f"0-{np.nanmax(p):.0f} dbar")

    CT_all, SA_all = np.array(CT_list), np.array(SA_list)
    ref_used = pd.DataFrame(used).reset_index(drop=True)
    n_ref = np.isfinite(CT_all).sum(axis=0)
    keep = n_ref >= REF_MIN_N
    with warnings.catch_warnings():                # levels without data -> NaN, no warning
        warnings.simplefilter("ignore", RuntimeWarning)
        CT_mean, SA_mean = np.nanmean(CT_all, axis=0), np.nanmean(SA_all, axis=0)
        CT_std, SA_std = np.nanstd(CT_all, axis=0), np.nanstd(SA_all, axis=0)

    P = ("pressure",)
    ds = xr.Dataset(
        {"CT": (P, CT_mean[keep], dict(long_name="mean Conservative Temperature", units="degC")),
         "SA": (P, SA_mean[keep], dict(long_name="mean Absolute Salinity", units="g kg-1")),
         "sig0": (P, gsw.sigma0(SA_mean, CT_mean)[keep],
                  dict(long_name="potential density anomaly (ref 0 dbar) of the mean CT, SA",
                       units="kg m-3")),
         "CT_std": (P, CT_std[keep], dict(long_name="std of CT across casts", units="degC")),
         "SA_std": (P, SA_std[keep], dict(long_name="std of SA across casts", units="g kg-1")),
         "n": (P, n_ref[keep].astype("int32"), dict(long_name="number of casts at this level")),
         "CT_cast": (("cast",) + P, CT_all[:, keep], dict(long_name="CT of each cast", units="degC")),
         "SA_cast": (("cast",) + P, SA_all[:, keep], dict(long_name="SA of each cast", units="g kg-1")),
         "station": ("cast", ref_used["station"].to_numpy("int32")),
         "lon": ("cast", ref_used["lon"].to_numpy(float), dict(units="degrees_east")),
         "lat": ("cast", ref_used["lat"].to_numpy(float), dict(units="degrees_north")),
         "time": ("cast", pd.to_datetime(ref_used["time"]).to_numpy("datetime64[ns]")),
         "outside_eddies": ("cast", ref_used["outside_eddies"].to_numpy(bool).astype("int8")),
         "file": ("cast", np.array(files, dtype=object))},
        coords={"pressure": ("pressure", p_ref[keep], dict(units="dbar", long_name="sea water pressure (bin centre)"))},
        attrs=dict(title="WHIRLS CTD background (reference) profile",
                   casts="all CTD casts" if use_all_casts else "CTD casts outside every TOEddies contour",
                   n_casts=len(ref_used), bin_dbar=REF_DP, min_casts_per_level=REF_MIN_N,
                   teos10="SA = gsw.SA_from_SP, CT = gsw.CT_from_t, downcasts only",
                   created_by="analysis/make_ctd_reference_profile.py",
                   created=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")))
    print(f"\nbackground from {len(ref_used)} casts, 0-{p_ref[keep].max():.0f} dbar "
          f"(levels with >= {REF_MIN_N} casts)")
    print(f"sigma0 range of the background: {float(ds.sig0.min()):.3f} .. {float(ds.sig0.max()):.3f}")
    return ds


def load_ctd_reference(path=OUT_ALL):
    """-> (ctd_ref, ref_used, CT_all, SA_all, p_ref), the same objects the notebooks used to build:
    ctd_ref = dict(p, CT, SA, sig0, CT_std, SA_std, n); ref_used = DataFrame(station, lon, lat, time)."""
    with xr.open_dataset(path) as ds:
        ds = ds.load()
    ctd_ref = dict(p=ds["pressure"].values, **{k: ds[k].values for k in
                                               ["CT", "SA", "sig0", "CT_std", "SA_std", "n"]})
    ref_used = pd.DataFrame({k: ds[k].values for k in ["station", "lon", "lat", "time", "outside_eddies"]})
    ref_used["outside_eddies"] = ref_used["outside_eddies"].astype(bool)
    return ctd_ref, ref_used, ds["CT_cast"].values, ds["SA_cast"].values, ds["pressure"].values


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outside", action="store_true",
                    help="only the casts outside every TOEddies contour")
    args = ap.parse_args()
    ds = make_reference(use_all_casts=not args.outside)
    out = OUT_OUTSIDE if args.outside else OUT_ALL
    ds.to_netcdf(out)
    print(f"saved -> {out}")
