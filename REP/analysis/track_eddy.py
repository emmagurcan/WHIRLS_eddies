import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import xarray as xr
import matplotlib.cm as cm
from matplotlib.colors import Normalize

import pyTOEddies_functions as f
import class_EDDIES as ed
from datetime import timedelta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def plot_mcells_modif(ax,cartmap,mat_eddies,type):
    """load eddy detection information from daily adt_yyyymmdd.mat
    and plot the eddy positions & contours
    """
    proj,crs,extent = cartmap
    xv,yv = f.box_square(extent[2],extent[-1],extent[0],extent[1])
    if type==-1:
        VAR,Name  = 'Cyclonic_Cell','Ncyclo'
        col  = 'b'
    elif type==1:
        VAR,Name  = 'Anticyclonic_Cell','Nanti'
        col='r'
    TCell = mat_eddies[VAR]
    Neddies=int(mat_eddies[Name][:][0][0])
    for n in range(0,Neddies):
        xcen = f.mod_shift(mat_eddies[TCell[2][n]][:])
        ycen = mat_eddies[TCell[3][n]][:]
        control = f.inpolygon(xcen,ycen, xv, yv)
        if control == True:
            xout = f.mod_shift(mat_eddies[TCell[4][n]][:])
            yout = mat_eddies[TCell[5][n]][:]
            xmax = f.mod_shift(mat_eddies[TCell[8][n]][:])
            ymax = mat_eddies[TCell[9][n]][:]
            ax.scatter(xcen,ycen,10,c='k',zorder=15,transform=crs)
            ax.plot(xmax,ymax,'-',transform=crs,color=col,linewidth=2,zorder=15)
            #ax.plot(xout,yout,'k--',transform=crs,linewidth=1,zorder=15)
          

# === Settings ===
dir_eddies = "/Users/emmagurcan/Documents/France/ENS/M1/stageM1/TOEddies/TOEddiesMSM40/Eddies/"
dir_bathy = "/Users/emmagurcan/Documents/France/ENS/M1/stageM1/data/bathymetry/etopo2.nc"


proj = ccrs.Mercator(central_longitude=-50, min_latitude=50, max_latitude=65)
crs = ccrs.PlateCarree()
pos_west,pos_east,pos_south,pos_north  = -43,-22, 56,62  #Reset the area
extent=[pos_west, pos_east, pos_south, pos_north]
# extent = [-45, -30, 56, 61]
cartmap = proj, crs, extent

bathymetry_data = xr.open_dataset(dir_bathy)
eddies = ed.Eddies(dir_eddies)

# Starting location
# MSM40 eddy 1 section 5
# current_x = -53.829277
# current_y = 59.4433842
# start_date = datetime(2014, 8, 18)
# end_date = datetime(2014, 8, 24)
# eddy_type = 1  # anticyclonic

# MSM40 eddy 2 section 5
# current_x = -52.24
# current_y = 59.34
# start_date = datetime(2014, 8, 19)
# end_date = datetime(2014, 8,31)
# eddy_type = 1

# MSM40 eddy 1 section 7
# current_x = -36.1
# current_y = 60.5
# start_date = datetime(2014, 8, 7)
# end_date = datetime(2014, 8, 24)
# eddy_type = -1

# MSM40 eddy 2 section 7
# current_x = -34.578679
# current_y = 58.901539
# start_date = datetime(2014, 8, 20)
# end_date = datetime(2014, 9, 2)
# eddy_type = -1

# # MSM40 eddy 2 section 7
# current_x = -30
# current_y = 58.17
# start_date = datetime(2014, 8, 15)
# end_date = datetime(2014, 8, 29)
# eddy_type = -1

# MSM40 eddy 2 section 7
current_x = -27
current_y = 58
start_date = datetime(2014, 6, 22)
end_date = datetime(2014, 9, 2)
eddy_type = 1
# Prepare to store centers and colors
eddy_centers = []
dates_used = []

# Generate date list
fx = current_x
fy = current_y


delta = timedelta(days=1)

date = start_date
t_cnt = 0
while date <= end_date:
    t_cnt += 1
    fname = f"adt_{date.strftime('%Y%m%d')}.mat"
    fpath = os.path.join(dir_eddies, fname)
    if not os.path.exists(fpath):
        date += delta
        continue

    print(f"Processing {fname}...")

    mat_eddies = h5py.File(fpath, "r")
    TCell = mat_eddies['Anticyclonic_Cell'] if eddy_type == 1 else mat_eddies['Cyclonic_Cell']
    Name = 'Nanti' if eddy_type == 1 else 'Ncyclo'
    Neddies = int(mat_eddies[Name][:][0][0])

    proj, crs, extent = cartmap
    xv, yv = f.box_square(extent[2], extent[-1], extent[0], extent[1])

    closest_dist = np.inf
    xcen_closest, ycen_closest = None, None

    for n in range(Neddies):
        xcen = f.mod_shift(mat_eddies[TCell[2][n]][:])
        ycen = mat_eddies[TCell[3][n]][:]
        control = f.inpolygon(xcen, ycen, xv, yv)
        if control:
            dist = np.hypot(xcen - current_x, ycen - current_y)
            if dist < closest_dist:
                closest_dist = dist
                xcen_closest, ycen_closest = xcen, ycen

    if xcen_closest is not None:
        eddy_centers.append((xcen_closest, ycen_closest))
        dates_used.append(date.strftime("%Y-%m-%d"))
        current_x, current_y = xcen_closest, ycen_closest

    date += delta

# === Plotting ===
fig = plt.figure(figsize=(8, 8), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=proj, zorder=3)
ax.add_feature(cfeature.LAND, linewidth=0.5, zorder=4, edgecolor='black', facecolor='gray')
ax.coastlines(resolution='10m', color='k', linestyle='-', alpha=1, zorder=5)
ax.set_extent(extent)

map_grid = ax.gridlines(linewidth=1, color='black', alpha=0.2, linestyle='--', draw_labels=True, zorder=6)
map_grid.top_labels = False
map_grid.left_labels = True
map_grid.right_labels = False
ft = 10
map_grid.xlabel_style = {'size': ft}
map_grid.ylabel_style = {'size': ft}
ax.set_extent(extent)
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='gray', zorder=3)
ax.coastlines(resolution='10m', color='k', linewidth=0.6, zorder=4)

# for filename in os.listdir(dir_eddies):
#     # Parse the date from the filename
#     date_str = filename.split("_")[1]
#     year = date_str[:4]
#     month = date_str[4:6]
#     day = date_str[6:8]
    
#     if day != '07' or month != '08' or year != '2014':
#         continue  # Skip all other dates

#     time_ed = datetime(int(year), int(month), int(day))
#     plot_mcells_modif(ax, cartmap, mat_eddies, -1)  # cyclonic

# Bathymetry
LON = bathymetry_data['lon']
LAT = bathymetry_data['lat']
depth = bathymetry_data['topo']
lon_subset = LON.sel(lon=slice(-62, -22))
lat_subset = LAT.sel(lat=slice(47, 65))
depth_subset = depth.sel(lon=slice(-62, -22), lat=slice(47, 65))


c = ax.contourf(lon_subset, lat_subset, depth_subset, np.arange(-4800, 500, 500),
                transform=crs, cmap='Blues_r', zorder=0, alpha=0.85)



# Plot eddy centers
colors = cm.autumn(np.linspace(0, 1, len(eddy_centers)))
for i, ((x, y), color) in enumerate(zip(eddy_centers, colors)):
    ax.scatter(x, y, color=color, s=12, transform=crs, zorder=10)
    # ax.text(x + 0.15, y, dates_used[i][5:], fontsize=5, transform=crs, color=color)
# ax.scatter(fx, fy, color='black', s=20, transform=crs, zorder=10)
title_str = f"Tracked Eddy Centers ({start_date.strftime('%b %d')} – {end_date.strftime('%b %d, %Y')})"
plt.title(title_str, fontsize=12)


# Create inset axis for horizontal colorbar just below the map
cax = inset_axes(ax,
                 width="80%",  # width relative to the parent axes
                 height="3%",  # height in %
                 loc='lower center',
                 bbox_to_anchor=(-0.105, -0.1, 1.2, 1),
                 bbox_transform=ax.transAxes,
                 borderpad=0)

cbar = fig.colorbar(c, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=10)
cbar.set_label('Depth [m]', fontsize=10, labelpad=5)
cbar.outline.set_linewidth(1)
# Optional: add a colorbar with a time scale
sm = plt.cm.ScalarMappable(cmap='autumn', norm=Normalize(vmin=0, vmax=len(dates_used)))
sm.set_array([])
cbar2 = fig.colorbar(sm, ax=ax, orientation='vertical', shrink=0.6)
cbar2.set_ticks(np.linspace(0, len(dates_used), 5))
cbar2.set_ticklabels([dates_used[int(i)] for i in np.linspace(0, len(dates_used) - 1, 5).astype(int)])
cbar2.set_label('Date progression', fontsize=9)

# plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.subplots_adjust(top=0.92, bottom=0.1, left=0.08, right=0.95)

plt.savefig("tracked_eddy_centers.png")
plt.close()

