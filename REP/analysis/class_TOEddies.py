import numpy as np
from matplotlib.path import Path

def inpolygon(xq, yq, xv, yv):
    shape = xq.shape
    xq = np.reshape(xq,-1)
    yq = np.reshape(yq,-1)
    xv = np.reshape(xv,-1)
    yv = np.reshape(yv,-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)


class Trajectories():
    "Tracking method"
    def __init__(self,name,mat,name_contour_max,data_max,name_contour_out,data_out,idname):
        "Define access to Trajectories and Contour data"
        self.name = name
        self.mat   = mat
        self.datamax = data_max
        self.name_contour_max = name_contour_max
        self.dataout = data_out
        self.name_contour_out = name_contour_out
        self.idname  = idname
        self.shapesmax= self.datamax[self.name_contour_max]
        self.shapesout= self.dataout[self.name_contour_out]
    def __str__(self):
        return 'Data for ({})'.format(self.name)
    def get_trajfield(self,number):
        '''Take chosen numbered field for all eddies
            -number (1-21) ex.id_traj = 0,date = 6 ...'''
        name,mat = self.name,self.mat
        traj  = mat[name]
        data =[];
        for column in traj[number,:]:
            value= mat[column][:]
            data.append(value)
        return data
    def get_trajfield_eddy(self,number,eddy):
        '''Take chosen field (number) for specific eddy'''
        name,mat = self.name,self.mat
        traj  = mat[name]
        data = mat[traj[number][eddy]][:]
        return data
    def get_trajxy(self,stps,eddies):
        '''Take the barycenters for eddies at given timestep (stps) along their whole trajectory'''
        name,mat = self.name,self.mat
        xbarys,ybarys,duration = [],[],[]
        for e in range(np.size(eddies,0)):
            stp = np.int(stps[e])
            eddy = np.int(eddies[e])
            dates  = self.get_trajfield_eddy(6,eddy)
            xibary = self.get_trajfield_eddy(4,eddy)#[stp]
            yibary = self.get_trajfield_eddy(5,eddy)#[stp]
            xbarys.append(xibary)
            ybarys.append(yibary)
            dura = dates[stp]-dates[0]
            duration.append(dura)
        return xbarys,ybarys,duration
    def find_idates(self,dates,py_date,xv,yv):
        '''Find the eddies that are in the box xv,yv in the given day'''
        N = np.size(dates)
        stps,eddies = [],[]
        for ii in range(N):
            estep = np.where(dates[ii]==py_date)
            if not np.size(estep)==0:
                xibary =  self.get_trajfield_eddy(4,ii)[estep]
                yibary =  self.get_trajfield_eddy(5,ii)[estep]
                logi=inpolygon(xibary,yibary, xv, yv)
                if logi ==True:
                    #eddies.append(ii)
                    #stps.append(estep[0])
                    eddies = np.append(eddies,ii)
                    stps = np.append(stps,estep[0])
        return(eddies,stps)
    def get_multishapes(self,stps,eddies,typecontour):
        '''I want give multi number of eddies and get the shapes in the day'''
        name,mat = self.name,self.mat
        idname = self.idname
        if typecontour =='max':
            data,field = self.datamax,self.name_contour_max
        elif typecontour =='out':
            data,field = self.dataout,self.name_contour_out
        shapes,id = data[field],np.array(data[idname])
        shapes_xmax,shapes_ymax=[],[]
        for e in range(np.size(eddies,0)):
            stp = np.int(stps[e])
            eddy = np.int(eddies[e])
            endex = mat[mat[name][1][eddy]][:][stp]
            [row,col] = np.where(id==endex);
            icol,irow =col[0],row[0]
            #shapes_xy.append(shapes[:,irow,icol])
            shapes_x,shapes_y = shapes[0,irow,icol], shapes[1,irow,icol]
            shapes_xmax.append(data[shapes_x][:])
            shapes_ymax.append(data[shapes_y][:])
        return shapes_xmax,shapes_ymax
