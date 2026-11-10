#!/usr/bin/env python3

"""This script creates the directory Slice_{i_slice}, containing a
time slice of instantaneous eddies. If Slice_{i_slice} already exists
then it is replaced.

Use this script if you have one file per date index, with the three
fields ssh, u and v together in this file, and you want to specify the
input files for all dates through a template file name.

Compared to inst_eddies.py, this script takes as arguments a template,
a first date index and a last date index, instead of a list of
files. So it is adequate when there is a date index in file
names. This script is useful if you do not want to list all the files
on the command line because the command line would be too long or if
some files are missing for some date indices.

A hyperslab is optionally selected with ncks.

The script reads the file inst_eddies_nml.txt, which must contain the
namelists `config_nml`, `input_ssh_nml` and
`read_grid_inst_eddies_nml`, in this order.

"""

import argparse
import sys
import glob
import os

import inst_eddies


def filename_generator(template, first_index, last_index):
    for i in range(first_index, last_index + 1):
        nc_file = template.format(i)

        # The template may contain a wild card as long as there is no
        # more than one file matching it. A missing file is ok and
        # will be dealt with in inst_eddies.py.

        path_list = glob.glob(nc_file)

        if len(path_list) == 1:
            nc_file = path_list[0]
        elif len(path_list) >= 2:
            sys.exit(
                "inst_eddies_L3.py:filename_generator: more than one "
                "file matching template"
            )

        yield nc_file


argparser = argparse.ArgumentParser()
argparser.add_argument(
    "--i_slice", type=int, help="slice number, integer", default=0
)
argparser.add_argument(
    "template",
    help="template of NetCDF file name (containing SSH and velocity at a "
    "single date), including path, with {} place holder in Python format "
    "specification mini-language (examples: {}, {:04}), and possibly shell "
    "wildcards in the rest of the path",
)
argparser.add_argument("first_index", type=int)
argparser.add_argument("last_index", type=int)
argparser.add_argument(
    "-b",
    "--bbox",
    nargs=4,
    type=float,
    metavar=("xmin", "xmax", "ymin", "ymax"),
)
args = argparser.parse_args()
inst_eddies.loop_inst_eddies(
    args.i_slice,
    args.first_index,
    filename_generator(args.template, args.first_index, args.last_index),
    args.bbox,
)
