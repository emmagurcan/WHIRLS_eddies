#!/bin/bash
# Run TOEddies end-to-end on one VarDyn-QG run (6-hourly, 0.1 deg).
#
#   mamba activate toeddies_pkgs
#   bash analysis/run_toeddies_vardyn.sh may      # or june / july
#
# Output goes to data/TOEDDIES/VARDYN_<run>/ :
#   Slice_0/            instantaneous eddies (shapefiles)
#   graph_C/, graph_AC/ overlap graph + trajectories for cyclones / anticyclones
#
# Time index ("date" in TOEddies) = number of 6-h steps since 2026-01-01 00h,
# so indices are comparable between the may/june/july runs:
#   date = 4*(day_of_year - 1) + hour/6     e.g. 2026-05-01 00h -> 480
set -e

RUN=${1:?usage: run_toeddies_vardyn.sh may|june|july}
ROOT=/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS
TOE=$ROOT/TOEddies
TRAJ=$TOE/Detection_eddies/Trajectories
RAW=$ROOT/data/VARDYN/VarDyn-QG-$RUN
CLEAN=$ROOT/data/VARDYN_toeddies/$RUN
OUT=$ROOT/data/TOEDDIES/VARDYN_$RUN

# Overlap parameters (in 6-h time steps and 0.1-deg grid points)
MAX_DELTA=4     # look up to 1 day ahead
DIST_LIM=10     # ~1 deg (~90-110 km) search radius for the extremum

# 1. Clean the files (lon/lat units, finite fill values)
python $ROOT/analysis/prep_vardyn_toeddies.py $RAW $CLEAN

# 2. Instantaneous eddies
mkdir -p $OUT && cd $OUT
cat > inst_eddies_nml.txt <<'EOF'
&config_nml min_amp = 0.001, max_radius = 50, 50, min_radius = 25./
&input_ssh_nml nc_time = 1, ssh_varname = "ssh", u_varname = "u", v_varname = "v"/
&read_grid_inst_eddies_nml periodic = f, uniform_lon_lat = t, stag_u = 0, 0, stag_v = 0, 0/
EOF
FIRST=$(ls $CLEAN/*.nc | head -1)
INIT=$(python - "$FIRST" <<'EOF'
import re, sys, datetime as dt
y, m, d, h = map(int, re.search(r"_y(\d{4})m(\d{2})d(\d{2})h(\d{2})", sys.argv[1]).groups())
print(4 * (dt.date(y, m, d) - dt.date(2026, 1, 1)).days + h // 6)
EOF
)
echo "init_date = $INIT"
$TOE/build/Inst_eddies/Create_slice/inst_eddies.py --init_date $INIT \
    $CLEAN/*.nc < inst_eddies_nml.txt

# 3. Overlap graph + trajectories, for each orientation
for ORIENT in C AC; do
    if [ $ORIENT = C ]; then CYC=t; NAME=Cyclones; else CYC=f; NAME=Anticyclones; fi
    cd $OUT
    printf '&config_graph_nml max_delta=%d, dist_lim=%d, cyclone=%s, min_inters_speed=0., min_inters_out=0.5 /\n&main_nml /\n' \
        $MAX_DELTA $DIST_LIM $CYC | $TOE/build/Overlap/eddy_graph.sh 1 .
    mkdir -p graph_$ORIENT
    mv edgelist.csv isolated_nodes.txt e_overestim.txt graph_$ORIENT/
    cd graph_$ORIENT
    $TRAJ/graph_properties.py ../ $NAME edgelist.csv output.gt
    $TRAJ/circumvent_nodes.py output.gt output1.gt
    $TRAJ/del_vert.py output1.gt output2.gt
    $TRAJ/delta_intra_seg.py output2.gt
    $TRAJ/stat_intra_seg.py delta_intra_seg.nc
    $TRAJ/cost_function.py output2.gt output3.gt
    $TRAJ/trajectories.py output3.gt output_final.gt
done

echo "Done: $OUT"