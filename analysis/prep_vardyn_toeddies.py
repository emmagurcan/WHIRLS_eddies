#!/usr/bin/env python3
"""Rewrite VarDyn-QG files in a form TOEddies' inst_eddies can read.

Why this is needed (VarDyn files as delivered):
  * lon/lat have no "units" attribute -> TOEddies (nf95_find_coord) cannot
    find the longitude/latitude coordinates and aborts.
  * _FillValue = NaN -> the "x == fill" test used to flag land never matches
    NaN, so land/missing points are not masked properly.
  * dimensions are called lon/lat; renaming to longitude/latitude also makes
    the optional --bbox (ncks hyperslab) of inst_eddies.py work.

Output: one file per time step, float32, with finite _FillValue and CF units.

Usage:
  python prep_vardyn_toeddies.py IN_DIR OUT_DIR [--ssh ssh] [--u ug] [--v vg]
"""
import argparse
import glob
import os
import re
from datetime import datetime

import netCDF4 as nc
import numpy as np

FILL = np.float32(1e20)

ap = argparse.ArgumentParser()
ap.add_argument("in_dir")
ap.add_argument("out_dir")
ap.add_argument("--ssh", default="ssh", help="SSH variable (ssh or sla)")
ap.add_argument("--u", default="ug", help="zonal velocity (ug or uc)")
ap.add_argument("--v", default="vg", help="meridional velocity (vg or vc)")
args = ap.parse_args()

os.makedirs(args.out_dir, exist_ok=True)
files = sorted(glob.glob(os.path.join(args.in_dir, "*.nc")))
if not files:
    raise SystemExit(f"no .nc files in {args.in_dir}")

pat = re.compile(r"_y(\d{4})m(\d{2})d(\d{2})h(\d{2})m(\d{2})\.nc$")

for f in files:
    m = pat.search(f)
    if not m:
        raise SystemExit(f"cannot parse date from {f}")
    y, mo, d, h, mi = map(int, m.groups())
    t = datetime(y, mo, d, h, mi)
    hours = (t - datetime(1950, 1, 1)).total_seconds() / 3600.0

    out = os.path.join(args.out_dir, os.path.basename(f))
    with nc.Dataset(f) as src, nc.Dataset(out, "w", format="NETCDF4") as dst:
        nlat, nlon = src.dimensions["lat"].size, src.dimensions["lon"].size
        dst.createDimension("time", 1)
        dst.createDimension("latitude", nlat)
        dst.createDimension("longitude", nlon)

        vt = dst.createVariable("time", "f8", ("time",))
        vt.units = "hours since 1950-01-01 00:00:00"
        vt.standard_name = "time"
        vt.calendar = "gregorian"
        vt[:] = hours

        vlon = dst.createVariable("longitude", "f4", ("longitude",))
        vlon.units = "degrees_east"
        vlon.standard_name = "longitude"
        vlon[:] = src["lon"][:]

        vlat = dst.createVariable("latitude", "f4", ("latitude",))
        vlat.units = "degrees_north"
        vlat.standard_name = "latitude"
        vlat[:] = src["lat"][:]

        for name_in, name_out, units in [
            (args.ssh, "ssh", "m"),
            (args.u, "u", "m s-1"),
            (args.v, "v", "m s-1"),
        ]:
            a = np.ma.filled(src[name_in][:].astype(np.float64), np.nan)
            a = np.where(np.isfinite(a), a, FILL).astype(np.float32)
            var = dst.createVariable(
                name_out, "f4", ("time", "latitude", "longitude"), fill_value=FILL
            )
            var.units = units
            var.source_variable = name_in
            var[:] = a

        dst.source = os.path.abspath(f)
        dst.history = "prep_vardyn_toeddies.py"

print(f"wrote {len(files)} files to {args.out_dir}")
