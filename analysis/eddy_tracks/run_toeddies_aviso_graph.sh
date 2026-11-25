#!/bin/bash
# Overlap graph + trajectories for the AVISO TOEddies run (all slices of data/TOEDDIES/AVISO_BIG).
#
#   bash analysis/run_toeddies_aviso_graph.sh C     # cyclones      -> AVISO_BIG/graph_C
#   bash analysis/run_toeddies_aviso_graph.sh AC    # anticyclones  -> AVISO_BIG/graph_AC
#
# The instantaneous eddies are not recomputed. AVISO_BIG holds one slice per period, in
# chronological order and without gaps (TOEddies requirement):
#   Slice_0  2026-02-01 -> 2026-03-31
#   Slice_1  2026-04-01 -> 2026-11-09
# To add an earlier period: detect it in a work folder with
#   inst_eddies_Aviso.py --i_slice 0 <template> <first> <last> < inst_eddies_nml.txt
# then rename the existing slices (Slice_k -> Slice_k+1), move the new one in as Slice_0,
# update n_slices.txt and rerun this script. A later period is simply Slice_<n>.
set -e

ORIENT=${1:?usage: run_toeddies_aviso_graph.sh C|AC}
export PATH=/Users/emmagurcan/miniforge3/envs/toeddies_pkgs/bin:$PATH
ROOT=/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS
TOE=$ROOT/TOEddies
TRAJ=$TOE/Detection_eddies/Trajectories
SHPC=$ROOT/data/TOEDDIES/AVISO_BIG

# Overlap parameters (daily steps, 1/8-deg grid points)
MAX_DELTA=5
DIST_LIM=32     # search radius for the extremum. Sweep on Feb-Mar 2026 (8..56): edges,
                # splits and merges saturate at 32 and are identical up to 56, so 32
                # truncates no overlapping pair (24 missed ~1% of splits/merges).

if [ $ORIENT = C ]; then CYC=t; NAME=Cyclones; else CYC=f; NAME=Anticyclones; fi
mkdir -p $SHPC/graph_$ORIENT && cd $SHPC/graph_$ORIENT
printf '&config_graph_nml max_delta=%d, dist_lim=%d, cyclone=%s, min_inters_speed=0., min_inters_out=0.5 /\n&main_nml /\n' \
    $MAX_DELTA $DIST_LIM $CYC | $TOE/build/Overlap/eddy_graph.sh 1 $SHPC
python $TRAJ/graph_properties.py $SHPC $NAME edgelist.csv output.gt
python $TRAJ/circumvent_nodes.py output.gt output1.gt
python $TRAJ/del_vert.py output1.gt output2.gt
python $TRAJ/delta_intra_seg.py output2.gt
python $TRAJ/stat_intra_seg.py delta_intra_seg.nc
python $TRAJ/cost_function.py output2.gt output3.gt
python $TRAJ/trajectories.py output3.gt output_final.gt
echo "Done: $SHPC/graph_$ORIENT"
