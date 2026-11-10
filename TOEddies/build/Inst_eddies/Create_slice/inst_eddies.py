#!/usr/bin/env python3

"""This script creates the directory Slice_{i_slice}, containing a
time slice of instantaneous eddies. If Slice_{i_slice} already exists
then it is replaced.

Use this script if you have one file per date, with the three fields
ssh, u and v together in this file, and you want to list the input
files for all dates (there can be any number of input files) on the
command line.

The script assumes that each input NetCDF file contains ssh, u and v
and reads ssh, u and v at the first date only of each input NetCDF
file. A hyperslab is optionally selected with ncks.

"""

import shutil
import os
import subprocess
import sys
import time
import csv

import nco
import shapefile
import f90nml

import create_slice


def process_1_file(
    bbox,
    nco_instance,
    nc_file,
    cropped_fname,
    inst_eddies_exe,
    writer,
    inst_eddies_nml,
    date,
):
    t0 = time.perf_counter()

    if bbox:
        xmin, xmax, ymin, ymax = bbox

        try:
            nco_instance.ncks(
                nc_file,
                output=cropped_fname,
                options=[
                    f"--dimension=longitude,{xmin},{xmax}",
                    f"--dimension=latitude,{ymin},{ymax}",
                ],
            )
        except nco.NCOException as err:
            # The error message coming from the binary NCO executable
            # is on standard output:
            print(err.stdout)
            raise

        nc_file = cropped_fname

    t1 = time.perf_counter()
    elapsed_NCO = t1 - t0
    t0 = t1
    subprocess.run(
        [inst_eddies_exe, nc_file],
        check=True,
        input=str(inst_eddies_nml) + "\n",
        text=True,
    )
    elapsed_Fortran = time.perf_counter() - t0
    writer.writerow([date, elapsed_NCO, elapsed_Fortran])


def loop_inst_eddies(i_slice, init_date, files, bbox):
    """`files` is an iterator of file paths. The first file must exist,
    others can be missing. bbox is a tuple of four real numbers.

    """

    print(
        "Enter the namelists `config_nml`, `input_ssh_nml` and "
        "`read_grid_inst_eddies_nml`, followed by CTRL-D:"
    )
    inst_eddies_nml = f90nml.read(sys.stdin)

    if "uniform_lon_lat" not in inst_eddies_nml["read_grid_inst_eddies_nml"]:
        inst_eddies_nml["read_grid_inst_eddies_nml"]["uniform_lon_lat"] = True
        # (same default value as in Fortran program)

    inst_eddies_exe = "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/Inst_eddies/inst_eddies"

    if not os.access(inst_eddies_exe, os.X_OK):
        sys.exit(
            "inst_eddies.loop_inst_eddies:"
            + inst_eddies_exe
            + " not found or not executable"
        )

    nco_instance = nco.Nco()
    create_slice.create_slice(
        i_slice,
        init_date,
        with_proj=not inst_eddies_nml["read_grid_inst_eddies_nml"][
            "uniform_lon_lat"
        ],
    )

    slice_dir = f"Slice_{i_slice}"
    perf_report = open(f"{slice_dir}/perf_report.csv", "w", newline="")
    writer = csv.writer(perf_report, lineterminator="\n")
    writer.writerow(["", "s", "s"])  # units
    writer.writerow(["date", "elapsed NCO", "elapsed Fortran"])  # long names

    # Short names (without blanks):
    writer.writerow(["date", "elapsed_NCO", "elapsed_Fortran"])

    cropped_fname = f"{slice_dir}/cropped.nc"

    # First file is out of the loop because of leave_trace:
    nc_file = next(files)
    print("inst_eddies.loop_inst_eddies: Date:", init_date)
    print("inst_eddies.loop_inst_eddies: Input file:", nc_file, flush=True)

    inst_eddies_nml = f90nml.Namelist(
        main_nml=f90nml.Namelist(),
        config_nml=inst_eddies_nml["config_nml"],
        input_ssh_nml=inst_eddies_nml["input_ssh_nml"],
        read_grid_inst_eddies_nml=inst_eddies_nml["read_grid_inst_eddies_nml"],
    )

    if i_slice == 0:
        inst_eddies_nml["main_nml"]["leave_trace"] = True
    else:
        inst_eddies_nml["main_nml"]["leave_trace"] = False

    inst_eddies_nml["main_nml"]["slice"] = i_slice
    process_1_file(
        bbox,
        nco_instance,
        nc_file,
        cropped_fname,
        inst_eddies_exe,
        writer,
        inst_eddies_nml,
        init_date,
    )
    date = init_date + 1

    if i_slice == 0:
        inst_eddies_nml["main_nml"]["leave_trace"] = False

    for nc_file in files:
        print("inst_eddies.loop_inst_eddies: Date:", date)
        print("inst_eddies.loop_inst_eddies: Input file:", nc_file, flush=True)

        if os.access(nc_file, os.F_OK):
            process_1_file(
                bbox,
                nco_instance,
                nc_file,
                cropped_fname,
                inst_eddies_exe,
                writer,
                inst_eddies_nml,
                date,
            )
        else:
            print("inst_eddies.loop_inst_eddies: Missing file:", nc_file)

            for orient in ["Cyclones", "Anticyclones"]:
                with open(
                    f"{slice_dir}/{orient}/ishape_last.txt", "a"
                ) as ishape_last, shapefile.Reader(
                    f"{slice_dir}/{orient}/" "extremum"
                ) as shpr:
                    n_shapes = len(shpr)
                    ishape_last.write(str(n_shapes - 1) + "\n")

            with open(f"{slice_dir}/date_range.txt", "w") as f_obj:
                f_obj.write(f"{init_date} {date + 1}\n")

        date += 1

    perf_report.close()

    if bbox:
        os.remove(cropped_fname)


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--i_slice", type=int, help="slice number, integer", default=0
    )
    argparser.add_argument(
        "--init_date",
        type=int,
        help="date index for the first file, integer",
        default=0,
    )
    argparser.add_argument(
        "file",
        nargs="+",
        help="NetCDF file containing SSH and velocity at a single date",
    )
    argparser.add_argument(
        "-b",
        "--bbox",
        nargs=4,
        type=float,
        metavar=("xmin", "xmax", "ymin", "ymax"),
    )
    args = argparser.parse_args()

    if not os.access(args.file[0], os.F_OK):
        sys.exit(args.file[0] + " not found")

    loop_inst_eddies(args.i_slice, args.init_date, iter(args.file), args.bbox)
