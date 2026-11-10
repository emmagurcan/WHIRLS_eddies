#!/usr/bin/env python3

"""This script creates the directory Slice_{args.i_slice}, containing
a time slice of instantaneous eddies. If Slice_{args.i_slice} already
exists then it is replaced.

The script assumes that for each variable, ssh, u and v, all the dates
to process are in a single NetCDF file. u and v can be in the same
file, and can be in the same file as ssh. So the script takes between
one and three files as command line arguments.

"""

import argparse
import shutil
import os
import subprocess
import sys
import time
import csv

import f90nml
import netCDF4

import create_slice
import find_coord


def process_1_date(
    inst_eddies_exe, NetCDF_files, inst_eddies_nml, writer, my_date
):
    t0 = time.perf_counter()
    subprocess.run(
        [inst_eddies_exe, *NetCDF_files],
        check=True,
        input=str(inst_eddies_nml) + "\n",
        text=True,
    )
    elapsed_Fortran = time.perf_counter() - t0
    writer.writerow([my_date, elapsed_Fortran])


argparser = argparse.ArgumentParser()
argparser.add_argument(
    "--i_slice", type=int, help="slice number, integer", default=0
)
argparser.add_argument("--init_date", type=int, default=0)
argparser.add_argument(
    "--init_nc_time",
    type=int,
    default=0,
    help="initial time index, 0-based, to read in the NetCDF files",
)
argparser.add_argument(
    "NetCDF_files",
    metavar="NetCDF-file",
    help="[1-3] files for SSH, u, v",
    nargs="+",
)
args = argparser.parse_args()

if len(args.NetCDF_files) >= 4:
    sys.exit("Specify between 1 and 3 NetCDF files.")

for fname in args.NetCDF_files:
    if not os.access(fname, os.F_OK):
        sys.exit(fname + " not found")

inst_eddies_exe = "/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/Inst_eddies/inst_eddies"

if not os.access(inst_eddies_exe, os.X_OK):
    sys.exit(
        "inst_eddies_INALT.py:"
        + inst_eddies_exe
        + " not found or not executable"
    )

with netCDF4.Dataset(args.NetCDF_files[0]) as SSH_dataset:
    time_dimension = find_coord.find_coord(SSH_dataset, "time")

    if time_dimension is None:
        n_dates = 1
    else:
        n_dates = SSH_dataset.dimensions[time_dimension].size

print(
    "Enter the namelists `config_nml`, `input_ssh_nml` and "
    "`read_grid_inst_eddies_nml`, followed with CTRL-D:"
)
# (`nc_time` should not be in the namelist `input_ssh_nml`, as it will
# be set below.)
inst_eddies_nml = f90nml.read(sys.stdin)

if "uniform_lon_lat" not in inst_eddies_nml["read_grid_inst_eddies_nml"]:
    inst_eddies_nml["read_grid_inst_eddies_nml"]["uniform_lon_lat"] = True
    # (same default value as in Fortran program)

create_slice.create_slice(
    args.i_slice,
    args.init_date,
    with_proj=not inst_eddies_nml["read_grid_inst_eddies_nml"][
        "uniform_lon_lat"
    ],
)
perf_report = open(f"Slice_{args.i_slice}/perf_report.csv", "w", newline="")
writer = csv.writer(perf_report, lineterminator="\n")
writer.writerow(["", "s"])  # units

# Short names (without blanks):
writer.writerow(["date", "elapsed_Fortran"])

inst_eddies_nml["input_ssh_nml"]["nc_time"] = args.init_nc_time + 1
# (index is 1-based in Fortran)

# First date is out of the loop because of leave_trace:
print("inst_eddies_INALT.py: Date:", args.init_date, flush=True)
inst_eddies_nml = f90nml.Namelist(
    main_nml=f90nml.Namelist(),
    config_nml=inst_eddies_nml["config_nml"],
    input_ssh_nml=inst_eddies_nml["input_ssh_nml"],
    read_grid_inst_eddies_nml=inst_eddies_nml["read_grid_inst_eddies_nml"],
)

if args.i_slice == 0:
    inst_eddies_nml["main_nml"]["leave_trace"] = True
else:
    inst_eddies_nml["main_nml"]["leave_trace"] = False

inst_eddies_nml["main_nml"]["slice"] = args.i_slice
process_1_date(
    inst_eddies_exe, args.NetCDF_files, inst_eddies_nml, writer, args.init_date
)

if args.i_slice == 0:
    inst_eddies_nml["main_nml"]["leave_trace"] = False

my_date = args.init_date

for inst_eddies_nml["input_ssh_nml"]["nc_time"] in range(
    args.init_nc_time + 2, n_dates + 1
):
    my_date += 1
    print("inst_eddies_INALT.py: Date:", my_date, flush=True)

    process_1_date(
        inst_eddies_exe, args.NetCDF_files, inst_eddies_nml, writer, my_date
    )

perf_report.close()
