RUNNING TOEddies for AVISO

0. mamba activate toeddies_pkgs
1. Build the project
in TOEddies
cmake -B build -S Detection_eddies -DCMAKE_BUILD_TYPE:STRING=Release -DFETCH:BOOL=ON -DBUILD_GDAL:BOOL=OFF -DBUILD_PROJ:BOOL=OFF
cmake --build build

2. Run on AVISO files
For AVISO daily files, create inst_eddies_nml.txt
for 1/8 deg
&config_nml min_amp = 0.001, max_radius = 40, 40, min_radius = 25./
&input_ssh_nml nc_time = 1, ssh_varname = "adt", u_varname = "ugos", v_varname = "vgos"/
&read_grid_inst_eddies_nml periodic = f, uniform_lon_lat = t stag_u = 0, 0, stag_v = 0, 0/

then run 

build/Inst_eddies/Create_slice/inst_eddies_Aviso.py /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/data/SSH/april_nov/2026-%m-%d_cmems_data.nc 2026-04-01 2026-11-09 < inst_eddies_nml.txt

3. Detection of overlapping contours
a. 
Run build/Overlap/eddy_graph.sh 1 . in TOEddies
for AVISO 1/8:
&CONFIG_GRAPH_NML
 MAX_DELTA=5          ,
 DIST_LIM=24         ,
 CYCLONE=T,
 MIN_INTERS_SPEED=  0.00000000    ,
 MIN_INTERS_OUT= 0.500000000    ,
 /
 For AC: Cyclone=F
 &CONFIG_GRAPH_NML
 MAX_DELTA=5          ,
 DIST_LIM=12         ,
 CYCLONE=T,
 MIN_INTERS_SPEED=  0.00000000    ,
 MIN_INTERS_OUT= 0.500000000    ,

&MAIN_NML /

mv edgelist.csv isolated_nodes.txt e_overestim.txt Cyclones 
Then run the same for anticyclones
mv edgelist.csv isolated_nodes.txt e_overestim.txt Anticyclones 

b. Run these files: 
graph_properties.py
circumvent_nodes.py
del_vert.py
delta_intra_seg.py
stat_intra_seg.py
cost_function.py
trajectories.py

cd to Cyclones and run graph_properties
/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/graph_properties.py ../ Cyclones edgelist.csv cyc_output.gt

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/circumvent_nodes.py cyc_output.gt cyc_output1.gt

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/del_vert.py cyc_output1.gt cyc_output2.gt

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/delta_intra_seg.py cyc_output2.gt

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/stat_intra_seg.py delta_intra_seg.nc

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/cost_function.py cyc_output2.gt cyc_output3.gt

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/trajectories.py cyc_output3.gt output_final.gt

to plot:
/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/Analysis/plot_traj.py expanded_traj.json ../
--annotate to annotate
--first_node xxx to show only 1 node

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/Analysis/draw_segments.py output_final.gt out.gv
xdot out.gv

For VarDyn
mkdir -p /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/data/TOEDDIES/VARDYN_july && cd /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/data/TOEDDIES/VARDYN_july

cat > inst_eddies_nml.txt <<'EOF'&config_nml min_amp = 0.001, max_radius = 50, 50, min_radius = 25./                                      
&input_ssh_nml nc_time = 1, ssh_varname = "ssh", u_varname = "u", v_varname = "v"/
&read_grid_inst_eddies_nml periodic = f, uniform_lon_lat = t, stag_u = 0, 0, stag_v = 0, 0/
EOF

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/Inst_eddies/Create_slice/inst_eddies.py --init_date 568 /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/data/VARDYN/VarDyn-QG-july/VarDyn-QG-july_y2026m*.nc < inst_eddies_nml.txt

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/build/Overlap/eddy_graph.sh 1 .

mkdir graph_C

mv edgelist.csv isolated_nodes.txt e_overestim.txt graph_C

cd graph_C                   

/Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/del_vert.py output1.gt output2.gt              
    /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/delta_intra_seg.py output2.gt
    /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/stat_intra_seg.py delta_intra_seg.nc
    /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/cost_function.py output2.gt output3.gt
    /Users/emmagurcan/Documents/France/ENS/Speich/WHIRLS/TOEddies/Detection_eddies/Trajectories/trajectories.py output3.gt output_final.gt
    
same for ac

For VarDyn: &CONFIG_GRAPH_NML MAX_DELTA=4, DIST_LIM=10 /
