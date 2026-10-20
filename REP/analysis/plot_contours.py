import matplotlib.pyplot as plt
import h5py
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os as os 
import numpy as np
import cmocean
import netCDF4 as nc
import fonctions as funcs

#For TOEddies
#import sys
#sys.path.insert(1, '\\TOEeddies\\TOEddiesMSM40')
import pyTOEddies_functions as f
import plot_toe_eddies as plot
import class_EDDIES as ed

from cartopy.feature import LAND
from cartopy.feature import ShapelyFeature
from cartopy import feature as cfeature
import xarray as xr

dir_eddies = "/Users/emmagurcan/Documents/France/ENS/M1/stageM1/TOEddies/TOEddiesMSM40/Eddies/"
dir_bathy = "/Users/emmagurcan/Documents/France/ENS/M1/stageM1/data/bathymetry/etopo2.nc"

bathymetry_data = xr.open_dataset(dir_bathy)
land_feature = cfeature.NaturalEarthFeature(
    category='physical',
    name='land',
    scale='10m',
    edgecolor='black',
    facecolor='gray'
)

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
           

ft=10

proj = ccrs.Mercator(central_longitude=-50,min_latitude=50,max_latitude=65)
crs = ccrs.PlateCarree()
pos_west,pos_east,pos_south,pos_north  = -62,-27, 47,62  #Reset the area
extent=[pos_west, pos_east, pos_south, pos_north]
cartmap = proj,crs,extent

# for filename in os.listdir(dir_eddies) :
#     matADT=filename
#     print(dir_eddies+matADT)
#     mat_eddies = h5py.File(dir_eddies+matADT,"r") #mat_eddies.keys()
#     time_ed = datetime(int(filename[4:8]), int(filename[8:10]), int(filename[10:12]))
#     date = time_ed.strftime("%m_%d")
    
#     fig = plt.figure(dpi=300)
#     ax0 = fig.add_subplot(1, 1, 1,projection=proj, zorder=3)
#     ax0.add_feature(land_feature, zorder=4)

#     # ax0.add_feature(cfeature.LAND,linewidth=0.5,zorder=4, edgecolor='black', facecolor='gray')#,alpha=0.1)
#     ax0.coastlines(resolution='10m', color='k', linestyle='-', alpha=1,zorder=5)
#     ax0.set_extent(extent)
#     map_grid=ax0.gridlines(linewidth=1, color='black', alpha=0.2, linestyle='--', draw_labels=True,zorder=6)
#     map_grid.top_labels = False
#     map_grid.left_labels = True
#     map_grid.right_labels = False
#     map_grid.xlabel_style = {'size': ft}
#     map_grid.ylabel_style = {'size': ft}

#     x, y,adt,u,v =eddies.get_data(time_ed)[:5]

#     #plot adt contours
#     #☺ax0.contour(x, y, adt, np.arange(-1,1,0.01),colors='k', transform=crs, linestyles='solid', linewidths=1, alpha=0.3)

#     pc = ax0.contourf(x, y, adt, np.arange(-1,0.71,0.01), cmap=cmocean.cm.dense, transform=crs, extent='both',  alpha=0.7)
#     cbar = fig.colorbar(pc)
#     cbar.ax.tick_params(labelsize=10)
#     cbar.set_label('ADT [m]',fontsize= 10 ,labelpad=10)
#     cbar.outline.set_linewidth(1)

#     #plot eddies contours
#     plot_mcells_modif(ax0,cartmap,mat_eddies,-1) #cyclonic
#     plot_mcells_modif(ax0,cartmap,mat_eddies,1) #anticyclonic
  
#     plt.savefig(f'{date}.png')
    
#     plt.close()
           
         
# for filename in os.listdir(dir_eddies) :
#     matADT=filename
#     mat_eddies = h5py.File(dir_eddies+matADT,"r") #mat_eddies.keys()
#     time_ed = datetime(int(filename[4:8]), int(filename[8:10]), int(filename[10:12]))

#     fig = plt.figure(dpi=300)

#     ax0 = fig.add_subplot(1, 1, 1,projection=proj, zorder=3)
#     ax0.set_title('01 juin 2014')
#     ax0.add_feature(cfeature.LAND,linewidth=0.5,zorder=4, edgecolor='black', facecolor='gray')#,alpha=0.1)
#     ax0.coastlines(resolution='10m', color='k', linestyle='-', alpha=1,zorder=5)
#     ax0.set_extent(extent)
#     map_grid=ax0.gridlines(linewidth=1, color='black', alpha=0.2, linestyle='--', draw_labels=True,zorder=6)
#     map_grid.top_labels = False
#     map_grid.left_labels = True
#     map_grid.right_labels = False
#     map_grid.xlabel_style = {'size': ft}
#     map_grid.ylabel_style = {'size': ft}

#     x, y,adt,u,v =eddies.get_data(time_ed)[:5]

#     #plot adt contours
#     #ax0.contour(x, y, adt, np.arange(0.1,0.8,0.01),colors='k', transform=crs, linestyles='solid', linewidths=1, alpha=0.3)

#     pc = ax0.contourf(x, y, adt,cmap=cmocean.cm.dense, transform=crs, extent='both',  alpha=0.7)
#     cbar = fig.colorbar(pc)
#     cbar.ax.tick_params(labelsize=10)
#     cbar.set_label('ADT [m]',fontsize= 10 ,labelpad=10)
#     cbar.outline.set_linewidth(1)

#     #ax0.scatter([-25],[14.5],color='yellow',transform=crs,zorder=14)

#     #plot eddies contours
#     plot_mcells_modif(ax0,cartmap,mat_eddies,-1) #cyclonic
#     plot_mcells_modif(ax0,cartmap,mat_eddies,1) #anticyclonic
           


# # Pour tout afficher

# ft=10

# #Set the projection
# proj = ccrs.Mercator(central_longitude=55,min_latitude=4,max_latitude=18)
# crs = ccrs.PlateCarree()
# pos_west,pos_east,pos_south,pos_north  = -62.5,-50,9,13  #Reset the area
# extent=[pos_west, pos_east, pos_south, pos_north]
# cartmap = proj,crs,extent

# eddies = ed.Eddies(dir_eddies)


# for filename in os.listdir(dir_eddies) :
#     matADT=filename
#     mat_eddies = h5py.File(dir_eddies+matADT,"r") #mat_eddies.keys()
#     time_ed = datetime(int(filename[4:8]), int(filename[8:10]), int(filename[10:12]))

#     fig = plt.figure()

#         #Salinity map
#     ax0 = fig.add_subplot(1, 1, 1,projection=proj, zorder=3)

#     ax0.add_feature(cfeature.LAND,linewidth=0.5,zorder=4, edgecolor='black', facecolor='darkseagreen')#,alpha=0.1)
#     ax0.coastlines(resolution='10m', color='k', linestyle='-', alpha=1,zorder=5)
#     ax0.set_extent(extent)
#     map_grid=ax0.gridlines(linewidth=1, color='black', alpha=0.2, linestyle='--', draw_labels=True,zorder=6)
#     map_grid.top_labels = False
#     map_grid.left_labels = True
#     map_grid.right_labels = False
#     map_grid.xlabel_style = {'size': ft}
#     map_grid.ylabel_style = {'size': ft}

#     x, y,adt,u,v =eddies.get_data(time_ed)[:5]
   
#     #plot adt contours
#     ax0.contour(x, y, adt, np.arange(0.1,0.8,0.01),colors='k', transform=crs, linestyles='solid', linewidths=1, alpha=0.3)
   
#     #plot eddies contours
#     plot_mcells_modif(ax0,cartmap,mat_eddies,-1) #cyclonic
#     plot_mcells_modif(ax0,cartmap,mat_eddies,1) #anticyclonic

import cmocean
import matplotlib.pyplot as plt
import h5py
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os as os 
import numpy as np
import cmocean
import xarray as xr


#For TOEddies
#import sys
#sys.path.insert(1, '\\TOEeddies\\TOEddiesMSM40')
import pyTOEddies_functions as f
import plot_toe_eddies as plot
import class_EDDIES as ed

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
            # ax.plot(xout,yout,'k--',transform=crs,linewidth=1,zorder=15)

ft=10

#Set the projection
proj = ccrs.Mercator(central_longitude=-50,min_latitude=50,max_latitude=65)
crs = ccrs.PlateCarree()
pos_west,pos_east,pos_south,pos_north  = -62,-20,47,62  #Reset the area
extent=[pos_west, pos_east, pos_south, pos_north]
cartmap = proj,crs,extent

# Load bathymetric data (example using NetCDF file)


# Extract relevant variables (e.g., longitude, latitude, depth)
lon = bathymetry_data['lon']
lat = bathymetry_data['lat']
depth = bathymetry_data['topo']

# Define region of interest
lon_min, lon_max, lat_min, lat_max = -62, -20, 47, 65

# Subset data to region of interest
lon_subset = lon.sel(lon=slice(lon_min, lon_max))
lat_subset = lat.sel(lat=slice(lat_min, lat_max))
depth_subset = depth.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))

eddies = ed.Eddies(dir_eddies)

filename="adt_20140815.mat"
mat_eddies = h5py.File(dir_eddies+filename,"r") #mat_eddies.keys()
time_ed = datetime(int(filename[4:8]), int(filename[8:10]), int(filename[10:12]))

fig = plt.figure(dpi=300)

ax0 = fig.add_subplot(1, 1, 1,projection=proj, zorder=3)
ax0.set_title('15 aout 2014')
ax0.add_feature(cfeature.LAND,linewidth=0.5,zorder=4, edgecolor='black', facecolor='gray')#,alpha=0.1)
ax0.coastlines(resolution='10m', color='k', linestyle='-', alpha=1,zorder=5)
ax0.set_extent(extent)
map_grid=ax0.gridlines(linewidth=1, color='black', alpha=0.2, linestyle='--', draw_labels=True,zorder=6)
map_grid.top_labels = False
map_grid.left_labels = True
map_grid.right_labels = False
map_grid.xlabel_style = {'size': ft}
map_grid.ylabel_style = {'size': ft}

# Plot bathymetric relief
c = ax0.contourf(lon_subset, lat_subset, depth_subset, np.arange(-4000,500,100), transform=ccrs.PlateCarree(), cmap='Blues_r',zorder=0,alpha=0.85)

x, y,adt,u,v =eddies.get_data(time_ed)[:5]

#plot adt contours
contour_lines = ax0.contour(x, y, adt, np.arange(-0.8, 0.1, 0.025), cmap='Greys', transform=crs, linestyles='dashed', linewidths=0.5)
ax0.clabel(contour_lines, inline=True, fmt='%1.1f', fontsize=7)
#pc = ax0.contour(x, y, adt,cmap=cmocean.cm.dense, transform=crs, extent='both',linestyles='dashed')
cbar = fig.colorbar(c)
cbar.ax.tick_params(labelsize=10)
cbar.set_label('Depth [m]',fontsize= 10 ,labelpad=10)
cbar.outline.set_linewidth(1)

ax0.scatter([-25],[14.5],color='yellow',transform=crs,zorder=14)

start_lon, start_lat = -60.0843, 59.8288  # Londres
end_lon, end_lat = -45.1687, 58.9023     # Paris

#Ajout de la ligne entre les deux points
ax0.plot([start_lon, end_lon], [start_lat, end_lat], color='black', linewidth=2, marker='o', transform=ccrs.Geodetic())

#plot eddies contours
plot_mcells_modif(ax0,cartmap,mat_eddies,-1) #cyclonic
plot_mcells_modif(ax0,cartmap,mat_eddies,1) #anticyclonic

#%%

import matplotlib.pyplot as plt
import h5py
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os as os 
import numpy as np
import cmocean
import class_EDDIES as ed
import xarray as xr
import pyTOEddies_functions as f
import scipy

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
            xmax = f.mod_shift(mat_eddies[TCell[8][n]][:])
            ymax = mat_eddies[TCell[9][n]][:]
            ax.scatter(xcen,ycen,10,c='k',zorder=15,transform=crs)
            ax.plot(xmax,ymax,'-',transform=crs,color=col,linewidth=2,zorder=15)
           

#Set the projection
proj = ccrs.Mercator(central_longitude=-50,min_latitude=50,max_latitude=65)
crs = ccrs.PlateCarree()
pos_west,pos_east,pos_south,pos_north = -62,-20,47,62  #Reset the area
extent=[pos_west, pos_east, pos_south, pos_north]
cartmap = proj,crs,extent

# Load Data

eddies = ed.Eddies(dir_eddies)

# Define the section (section 2 can't be plotted => no uCTD measures)
rep = '/Users/emmagurcan/Documents/France/ENS/M1/stageM1'
uCTD = scipy.io.loadmat(rep+'/Data/Data_MSM40/UCTD_dataproc.mat')
adcp = scipy.io.loadmat(rep +'/Data/Data_MSM40/adcp_150kHz/TL_GEOMAR_150_STJOHNS_BREST_000_000000_6_hc.mat')
section_num = 7
# sub_section = 1
start_index, end_index, len_section, start_index_adcp, end_index_adcp = funcs.section(section_num)
lon = np.squeeze(uCTD['uctdlon'])[start_index:end_index+1]
lat = np.squeeze(uCTD['uctdlat'])[start_index:end_index+1]

for filename in os.listdir(dir_eddies) :
    matADT=filename
    mat_eddies = h5py.File(dir_eddies+matADT,"r") #mat_eddies.keys()
    time_ed = datetime(int(filename[4:8]), int(filename[8:10]), int(filename[10:12]))
    date = time_ed.strftime("%m_%d")
    
    fig = plt.figure(dpi=300)
    ax0 = fig.add_subplot(1, 1, 1,projection=proj, zorder=3)
    ax0.add_feature(cfeature.LAND,linewidth=0.5,zorder=4, edgecolor='black', facecolor='gray')#,alpha=0.1)
    ax0.coastlines(resolution='10m', color='k', linestyle='-', alpha=1,zorder=5)
    ax0.set_extent(extent)

    date_str = filename.split("_")[1]
    year = date_str[:4]
    month = date_str[4:6]
    day = date_str[6:8]

    ax0.set_title(f"{day}-{month}-{year}")
    map_grid=ax0.gridlines(linewidth=1, color='black', alpha=0.2, linestyle='--', draw_labels=True,zorder=6)
    map_grid.top_labels = False
    map_grid.left_labels = True
    map_grid.right_labels = False
    ft=10
    map_grid.xlabel_style = {'size': ft}
    map_grid.ylabel_style = {'size': ft}
    
    # Bathymetry
    LON = bathymetry_data['lon']
    LAT = bathymetry_data['lat']
    depth = bathymetry_data['topo']
    lon_min, lon_max, lat_min, lat_max = -62, -20, 47, 65
    lon_subset = LON.sel(lon=slice(lon_min, lon_max))
    lat_subset = LAT.sel(lat=slice(lat_min, lat_max))
    depth_subset = depth.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
    c = ax0.contourf(lon_subset, lat_subset, depth_subset, np.arange(-4000,500,500), transform=ccrs.PlateCarree(), cmap='Blues_r',zorder=0,alpha=0.85)
    cbar = fig.colorbar(c, ax=ax0, shrink=0.7)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('Depth [m]',fontsize= 10 ,labelpad=10)
    cbar.outline.set_linewidth(1)

    x, y,adt,u,v =eddies.get_data(time_ed)[:5]
    
    #plot adt contours
    contour_lines = ax0.contour(x, y, adt, np.arange(-0.8, 0.1, 0.02), cmap='Greys', transform=crs, linestyles='dashed', linewidths=0.7)
    #ax0.clabel(contour_lines, inline=True, fmt='%1.1f', fontsize=7)
    
    #plot eddies contours
    plot_mcells_modif(ax0,cartmap,mat_eddies,-1) #cyclonic
    plot_mcells_modif(ax0,cartmap,mat_eddies,1) #anticyclonic
    
    # Plot the different measurements and moorings
    plt.scatter(lon, lat, marker='o', color='y', s=5, label='uCTD', transform=ccrs.Geodetic(),zorder=20)
    
    plt.savefig(f'{date}.png')
    
    plt.close()
