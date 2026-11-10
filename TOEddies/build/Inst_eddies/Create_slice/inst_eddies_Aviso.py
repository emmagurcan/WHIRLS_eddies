#!/usr/bin/env python3

"""This script creates the directory Slice_{i_slice}, containing a
time slice of instantaneous eddies. If Slice_{i_slice} already exists
then it is replaced.

Use this script if you have one file per calendar date, with the three
fields ssh, u and v together in this file, and you want to specify the
input files for all dates through a template file name.

Compared to inst_eddies.py, this script takes as arguments a template,
a first calendar date and a last calendar date, included, instead of a
list of files. So it is adequate when there is a calendar date pattern
in file names.

A hyperslab is optionally selected with ncks.

The script reads the file inst_eddies_nml.txt.

"""

import argparse
import datetime
import sys
import glob
import os

from dateutil import parser

import inst_eddies


def filename_generator(template, first_date, last_date):
    my_date = first_date

    while my_date <= last_date:
        nc_file = my_date.strftime(template)

        # The template may contain a wild card as long as there is no
        # more than one file matching it. A missing file is ok and
        # will be dealt with in inst_eddies.py.

        path_list = glob.glob(nc_file)

        if len(path_list) == 1:
            nc_file = path_list[0]
        elif len(path_list) >= 2:
            sys.exit(
                "inst_eddies_Aviso.py:filename_generator: more than one "
                "file matching template"
            )

        yield nc_file
        my_date += datetime.timedelta(1)


argparser = argparse.ArgumentParser()
argparser.add_argument(
    "template",
    help="template of NetCDF file name (containing SSH and velocity at a "
    "single date), including path, with date.strftime format codes (and "
    "possibly shell wildcards)",
)
argparser.add_argument(
    "--i_slice", type=int, help="slice number, integer", default=0
)
argparser.add_argument("first_date", help="%%Y-%%m-%%d")
argparser.add_argument("last_date", help="%%Y-%%m-%%d, included")
argparser.add_argument(
    "-b",
    "--bbox",
    nargs=4,
    type=float,
    metavar=("xmin", "xmax", "ymin", "ymax"),
)
args = argparser.parse_args()
first_date = parser.parse(args.first_date).date()
last_date = parser.parse(args.last_date).date()
d = (first_date - datetime.date(1950, 1, 1)).days

# We choose to use d as the first date index. We could choose anything.
inst_eddies.loop_inst_eddies(
    args.i_slice,
    d,
    filename_generator(args.template, first_date, last_date),
    args.bbox,
)
