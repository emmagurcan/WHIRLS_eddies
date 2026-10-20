import numpy as np
import h5py

class Eddies():
    "AVISO"
    def __init__(self,dir_mat):
        self.dir_mat = dir_mat
    def get_data(self,time):
        "load AVISO information from daily adt_yyyymmdd.mat"
        #print(aviso.keys())
        matADT='adt_%s.mat'  %(time.strftime("%Y%m%d"))
        mat_eddies = h5py.File(self.dir_mat+matADT,"r")
        date = time.strftime("%d-%b-%Y")
        lon,lat  =  mat_eddies['X'][:],mat_eddies['Y'][:]
        lons, lats = np.meshgrid(lon,lat)
        adt = mat_eddies['ADT'][:] #
        u   = mat_eddies['U'][:]
        v   = mat_eddies['V'][:]
        data = lons, lats,adt,u,v,date
        return data