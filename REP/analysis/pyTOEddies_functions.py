#export PATH=/opt/anaconda3/bin/:$PATH
#export PATH=/data/guez/miniconda3/envs/cartopy/bin:$PATH


import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import gsw
from scipy.signal import argrelextrema
from scipy import integrate
from datetime import timedelta,datetime,date


def matlab2datetime(matlab_datenum):
    "convert MATLAB time to python time"
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return day + dayfrac


def mod_date(vector):
    time_vec=[];
    for p in range(np.size(vector)):
        if vector[p]>0:
            time_vec=np.append(time_vec,matlab2datetime(np.int(vector[p])))
        else:
            time_vec=np.append(time_vec,np.nan)
    return(time_vec)

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + n*timedelta(days=1)

def mod_shift(lon):
    "Shift longitude from [0,360] to [-180,180]"
    lon_shift = np.mod(lon+180,360.)-180
    return lon_shift
   
def mod_shift2(lon):
    "Shift longitude from [0,360] to [-180,180]"
    lon_shift = str(np.mod(lon+180,360.)-180)
    return lon_shift

def find_lifetime(dates,life_time):
    N = np.size(dates)
    id_no=[]
    id_ok=[]
    for ii in range(N):
        if np.size(dates[ii])<life_time:
            id_no.append(ii)
        else:
            id_ok.append(ii)
    #print(id_list)
    dates = np.delete(dates, id_no, axis=0)
    return(dates,id_ok,id_no)


#% TOEddies Trajectories
#%------------------------------------------------------%
#% "Filter trajectories based on their Lifetime"
#%------------------------------------------------------%
#%

def get_trajfilt(traj,weeks,days_in_week, length_date):
    life=np.size(length_date, 0)
    cut_off = days_in_week*weeks
    filt = np.where(life>=cut_off)[0]
    traj_filt = traj[filt,:]
    return traj_filt


def fix_lons(lons):
    """
    Fix the given longitudes into the range ``[-180, 180]``.

    """
    lons = np.array(lons, copy=False, ndmin=1)
    fixed_lons = ((lons + 180) % 360) - 180
    # Make the positive 180s positive again.
    fixed_lons[(fixed_lons == -180) & (lons > 0)] *= -1
    return fixed_lons


def cm_to_inch(value):
    return value/2.54

def inpolygon(xq, yq, xv, yv):
    "MATLAB inpolygon equivelant in python"
    from matplotlib.path import Path
    shape = xq.shape
    xq = np.reshape(xq,-1)
    yq = np.reshape(yq,-1)
    xv = np.reshape(xv,-1)
    yv = np.reshape(yv,-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)
   
   
def inpolygon_new(xq, yq,p):
    "MATLAB inpolygon equivelant in python"
    shape = xq.shape
    xq = np.reshape(xq,-1)
    yq = np.reshape(yq,-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    return p.contains_points(q).reshape(shape)

def box_square(llcrnrlat,urcrnrlat,llcrnrlon,urcrnrlon):
    yd,yu,xr,xl = llcrnrlat,urcrnrlat,llcrnrlon,urcrnrlon
    poly= [[xl,yd],
           [xr,yd],
           [xr,yu],
           [xl,yu],
           [xl,yd]]
    xv = [point[0] for point in poly]
    yv = [point[1] for point in poly]
    return(xv,yv)

def distance(s_lat, s_lng, e_lat, e_lng):
   R = 6373.0
   s_lat = s_lat*np.pi/180.0                     
   s_lng = np.deg2rad(s_lng)    
   e_lat = np.deg2rad(e_lat)                      
   e_lng = np.deg2rad(e_lng) 
   d = np.sin((e_lat - s_lat)/2)**2 + np.cos(s_lat)*np.cos(e_lat) * np.sin((e_lng - s_lng)/2)**2
   return 2 * R * np.arcsin(np.sqrt(d))





def get_MLD(D,Z,delta_rho):
    """Compute MLD with threshold method on density
    (delta_rho = 0.03 kgm^-3)
    """
    depth = np.abs(Z)
    refdep=10;
    I_10=np.max(np.where(depth<=refdep))
    mld=None;
    ref_pdens=D[I_10];    #density at refdep
    drho_profile = np.abs(D-ref_pdens)-delta_rho;
    purged_drho_profile  = drho_profile[np.logical_not(np.isnan(drho_profile))]
    if np.size(purged_drho_profile) == 0:
        mld = np.NaN
    else :
        try:
            ttt=np.min(np.where(drho_profile>0))
        except ValueError:  #raised if array is empty.
            ttt = None
            print('ttt =', ttt)
            print('size of ttt =', np.size(ttt))
        if ttt==None:
            mld = np.NaN
        else :
            mld = depth[np.max(np.where(drho_profile<=0))]
            mld = mld+refdep;  
    return(mld)