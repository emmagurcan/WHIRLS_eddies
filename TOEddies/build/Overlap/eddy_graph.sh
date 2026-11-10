#!/bin/bash

# This script creates the edge list for a given orientation of
# eddies. The script takes an SHPC as argument.

if (($# <= 1))
then
    echo "Usage: eddy_graph.sh N_PROC SHPC-directory"
    exit 1
fi

set -e
##set -x
echo "mpiexec: $mpiexec"
${mpiexec:-mpiexec} -n $1 /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/Overlap/eddy_graph $2

rm -f isolated_nodes.txt edgelist_unsorted.csv
rank=0

while ((rank < $1))
do
    cat edgelist_$rank.csv >>edgelist_unsorted.csv
    cat isolated_nodes_$rank.txt >>isolated_nodes.txt
    rm edgelist_$rank.csv isolated_nodes_$rank.txt
    ((rank += 1))
done

# Sort the edge list for easier comparison between runs:
export LC_NUMERIC=C # for the sort command
sort --key=1,1 --key=2,2 --numeric-sort --output=edgelist.csv \
     edgelist_unsorted.csv
rm edgelist_unsorted.csv

echo "Created files \"edgelist.csv\", \"isolated_nodes.txt\" and " \
     "\"e_overestim.txt\"."
