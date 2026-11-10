#!/usr/bin/env python3

import os
import pathlib
import shutil

import shapefile


def create_slice(i_slice, init_date, with_proj):
    if os.access(f"Slice_{i_slice}", os.F_OK):
        shutil.rmtree(f"Slice_{i_slice}")

    for orient in ["Cyclones", "Anticyclones"]:
        os.makedirs(f"Slice_{i_slice}/{orient}")

        with shapefile.Writer(
            f"Slice_{i_slice}/{orient}/extremum", shapefile.POINT
        ) as w:
            w.field("ssh", "N", 13, 6)
            w.field("date", "N", 5)
            w.field("eddy_index", "N", 5)
            w.field("speed", "N", 13, 6)

        if with_proj:
            with shapefile.Writer(
                f"Slice_{i_slice}/{orient}/extr_proj", shapefile.POINT
            ) as w:
                w.field("date", "N", 5)
                w.field("eddy_index", "N", 5)

        for fname in ["outermost_contour", "max_speed_contour"]:
            with shapefile.Writer(
                f"Slice_{i_slice}/{orient}/{fname}", shapefile.POLYGON
            ) as w:
                w.field("r_eq_area", "N", 10, 4)
                w.field("ssh", "N", 13, 6)
                w.field("date", "N", 5)
                w.field("eddy_index", "N", 5)

        pathlib.Path(f"Slice_{i_slice}/{orient}/ishape_last.txt").touch()

    pathlib.Path(f"Slice_{i_slice}/timings.txt").touch()

    with open(f"Slice_{i_slice}/date_range.txt", "w") as f_obj:
        f_obj.write(f"{init_date} {init_date}\n")


if __name__ == "__main__":
    import sys

    import f90nml

    if len(sys.argv) != 2:
        sys.exit("Required argument: initial-date")

    if os.access("inst_eddies_nml.txt", os.R_OK):
        print("create_slice.py: Will use inst_eddies_nml.txt")
    else:
        sys.exit(
            'create_slice.py: "inst_eddies_nml.txt" not found -- Aborting'
        )

    inst_eddies_nml = f90nml.read("inst_eddies_nml.txt")

    try:
        i_slice = inst_eddies_nml["main_nml"]["slice"]
    except KeyError:
        i_slice = 0  # same default value as in Fortran program

    try:
        with_proj = not inst_eddies_nml["read_grid_inst_eddies_nml"][
            "uniform_lon_lat"
        ]
    except KeyError:
        with_proj = False  # same default value as in Fortran program

    create_slice(i_slice=i_slice, init_date=sys.argv[1], with_proj=with_proj)
