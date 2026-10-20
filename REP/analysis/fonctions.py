import numpy as np
import gsw
import glob
import pandas as pd
from scipy import interpolate
from scipy.interpolate import griddata
from scipy import signal


# Fonction pour charger des données argo et calculer un profil moyen régional
# Le mois choisi ici pour la moyenne climatologique est aoüt
def argo(name, mnth="08"):
   
    # Définition des listes
    file_list = []
    Temperature = []
    Salinite = []
    Lon = []
    Lat = []
    P_argo = []
   
    # Chargement des données
    for filename in glob.glob(name + '\\*.csv'):
        file_list.append(filename)

    for i in file_list:
        data = pd.read_csv(i, sep=",")

        try:
            if "TEMP (degree_Celsius)" in data.columns and "PSAL (psu)" in data.columns and "LONGITUDE (degree_east)" in data.columns and "LATITUDE (degree_north)" in data.columns:
                Temperature = np.vstack(
                    (Temperature, data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7] == mnth, "TEMP (degree_Celsius)"].values[:, None]))
                Salinite = np.vstack(
                    (Salinite, data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7] == mnth, "PSAL (psu)"].values[:, None]))
                Lon = np.vstack((Lon, data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                                == mnth, "LONGITUDE (degree_east)"].values[:, None]))
                Lat = np.vstack((Lat, data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                                == mnth, "LATITUDE (degree_north)"].values[:, None]))
                P_argo = np.vstack(
                    (P_argo, data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7] == mnth, "PRES (decibar)"].values[:, None]))
        except:
            if "TEMP (degree_Celsius)" in data.columns and "PSAL (psu)" in data.columns and "LONGITUDE (degree_east)" in data.columns and "LATITUDE (degree_north)" in data.columns:
                Temperature = data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                                       == mnth, "TEMP (degree_Celsius)"].values[:, None]
                Salinite = data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                                    == mnth, "PSAL (psu)"].values[:, None]
                Lon = data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                               == mnth, "LONGITUDE (degree_east)"].values[:, None]
                Lat = data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                               == mnth, "LATITUDE (degree_north)"].values[:, None]
                P_argo = data.loc[data["DATE (YYYY-MM-DDTHH:MI:SSZ)"].str[5:7]
                                  == mnth, "PRES (decibar)"].values[:, None]

    # Calcul de la salinité absolue et de la température potentielle
    SA = gsw.SA_from_SP(Salinite, P_argo, Lon, Lat)
    Temperature = gsw.CT_from_t(SA, Temperature, P_argo)
    Sig_argo = gsw.sigma0(SA, Temperature)
   
    # Calcul des profils moyens régionaux

    # Définition des différents profils Argo
    ind_P_max = np.arange(int(np.nanmin(P_argo)), int(np.nanmax(P_argo)), 1)
    ind_stop_profile = [0]
    for i in range(P_argo.shape[0]-1):
        if abs(P_argo[i+1]-P_argo[i]) > 1500:
            ind_stop_profile.append(i)

    del ind_stop_profile[1:3]

    # Définition des matrices de température, salinité, profondeur et densité
    MAT_temp = np.zeros((len(ind_stop_profile), len(ind_P_max)))
    MAT_S = np.zeros((len(ind_stop_profile), len(ind_P_max)))
    MAT_temp[:] = np.nan
    MAT_S[:] = np.nan
    Depth = np.zeros((len(ind_stop_profile), len(ind_P_max)))
    Depth[:] = np.nan
    MAT_Sig = np.zeros((len(ind_stop_profile), len(ind_P_max)))
    MAT_Sig[:] = np.nan

    # Remplissage des matrices
    for i in range(1, len(ind_stop_profile)):
        ii = ind_stop_profile[i-1]
        iii = ind_stop_profile[i]

        Temp = Temperature[ii:iii, 0]
        S = SA[ii:iii, 0]
        P = P_argo[ii:iii, 0]
        lon = Lon[ii:iii, 0][0]
        lat = Lat[ii:iii, 0][0]
        Sig = Sig_argo[ii:iii, 0]

        arg_sort = np.argsort(P)
        P = P[arg_sort]
        Temp = Temp[arg_sort]
        S = S[arg_sort]
        Z = gsw.z_from_p(P, lat)
        Sig = Sig[arg_sort]

        ind_Z = np.arange(int(np.nanmax(Z))-1, int(np.nanmin(Z))+1, -1)

        f = interpolate.interp1d(Z, Temp, 'linear')
        MAT_temp[i, -ind_Z] = f(ind_Z)
        g = interpolate.interp1d(Z, S, 'linear')
        MAT_S[i, -ind_Z] = g(ind_Z)
        h = interpolate.interp1d(Z, Sig, 'linear')
        MAT_Sig[i, -ind_Z] = h(ind_Z)
        Depth[i, -ind_Z] = ind_Z

    # Calcul profondeur, température, salinité, et densité moyenne
    Z_bar = np.nanmean(Depth, axis=0)
    T_bar = np.nanmean(MAT_temp, axis=0)
    S_bar = np.nanmean(MAT_S, axis=0)
    Sig_bar = np.nanmean(MAT_Sig, axis=0)
   
    # Ecart-types température, salinité et densité
    std_T = np.nanstd(MAT_temp, axis=0)
    std_S = np.nanstd(MAT_S, axis=0)
    std_Sig = np.nanstd(MAT_Sig, axis=0)

    return T_bar,S_bar,Sig_bar,std_T,std_S,std_Sig



# Fonction qui calcule les anomalies (ici cas température qui peut être remplacé par salinité)
def anom(PT,T_moy,PD,Sig_moy,ind_rho,ind_X,ind_Z,depth):
   
    # TEMPERATURE ANOMALY SUR LES ISOPYCNES PUIS SUR LES NIVEAUX EN Z

    # Interpolation de la moyenne sur les isopycnes
    T_moy_interp = np.zeros(len(ind_rho))
    check = np.where(np.isnan(T_moy)==0)[0]
    ind = np.arange(np.where(ind_rho>=np.nanmin(Sig_moy[check]))[0][0],np.where(ind_rho<=np.nanmax(Sig_moy[check]))[0][-1])
    f = interpolate.interp1d(Sig_moy[check],T_moy[check],'linear')
    T_moy_interp[ind] = f(ind_rho[ind])

    # Variables interpolée sur les isopycnes
    len_section=len(ind_X)
    T_interp = np.zeros((len(ind_rho),len_section))
    T_interp[:] = np.nan
    Anom_T_rho = np.zeros((len(ind_rho),len_section))
    Anom_T_rho[:] = np.nan
    depT = np.zeros((len(ind_rho),len_section))
    depT[:] = np.nan
    Anom_T = np.zeros((len(ind_Z),len_section)) #sur les niveaux en Z
    Anom_T[:] = np.nan

    for i in ind_X:
        # Interpolation de la température sur les isopycnes
        check = np.where(np.isnan(PT[:,i])==0)[0]
        if len(check)>1 :
            ind2 = np.arange(np.where(ind_rho>=np.nanmin(PD[check,i]))[0][0],np.where(ind_rho<=np.nanmax(PD[check,i]))[0][-1])
            f = interpolate.interp1d(PD[check,i],PT[check,i],'linear')
            T_interp[ind2,i] = f(ind_rho[ind2])
            # Calcul de l'anomalie
            for index in ind2:
                if T_moy_interp[index] != 0:
                    Anom_T_rho[index,i] = T_interp[index,i] - T_moy_interp[index]
            # Interpolation de la profondeur sur les isopycnes
            g = interpolate.interp1d(PD[check,i],depth[check,i],'linear')
            depT[ind2,i] = g(ind_rho[ind2])
           
        # Interpolation des anomalies isopycnale sur l'interpolation précédente (anomalies en z)
        check = np.where(np.isnan(depT[:,i])==0)[0]
        if len(check)>1 :
            ind3 = np.arange(np.where(ind_Z<=np.nanmax(depT[check,i]))[0][0],np.where(ind_Z>=np.nanmin(depT[check,i]))[0][-1])
            m = interpolate.interp1d(depT[check,i],Anom_T_rho[check,i],'linear')
            Anom_T[ind3,i] = m(ind_Z[ind3])

    return Anom_T




# Fonction qui sélectionne la section chosie
def section(section):
    # Sélection de la section
    if (section==1):
        start_index = 0
        end_index = 21
        start_index_adcp = 0
        end_index_adcp = 1898
    elif (section==2):
        print("Pas de mesures uCTD")
        return 0
    elif (section==3):
        start_index = 24
        end_index = 44
        start_index_adcp = 8486
        end_index_adcp = 11324
    elif (section==4):
        start_index = 45
        end_index = 65
        start_index_adcp = 11325
        end_index_adcp = 13091
    elif (section==5):
        start_index = 66
        end_index = 105
        start_index_adcp = 13092
        end_index_adcp = 15924
    elif (section==6):
        start_index = 106
        end_index = 119
        start_index_adcp = 15925
        end_index_adcp = 18186
    elif (section==7):
        start_index = 120
        end_index = 165
        start_index_adcp = 18187
        end_index_adcp = 25525
   
    # Taille de la section
    len_section = end_index - start_index + 1
   
    return start_index,end_index,len_section,start_index_adcp,end_index_adcp

def interpolation2D(values, x, z, max_x, delta_x, max_z, delta_z):
    """
    Performs 2D interpolation on a given dataset, filling NaN values using 
    nearest neighbor interpolation and then interpolating onto a new regular grid.

    Parameters:
    ----------
    values : np.ndarray
        2D array of data values to be interpolated (can contain NaNs).
    x : np.ndarray
        1D array of x-coordinates corresponding to the columns of `values`.
    z : np.ndarray
        1D array of z-coordinates corresponding to the rows of `values`.
    max_x : float
        Maximum x value for the new grid.
    delta_x : float
        Grid spacing in the x-direction for the new interpolated grid.
    max_z : float
        Maximum z value for the new grid.
    delta_z : float
        Grid spacing in the z-direction for the new interpolated grid.

    Returns:
    -------
    new_values_reshaped : np.ndarray
        Interpolated 2D array on the new grid.
    new_distance : np.ndarray
        1D array of new x-coordinates for the interpolated grid.
    new_depth : np.ndarray
        1D array of new z-coordinates for the interpolated grid.
    """
    
    # Find indices of NaN values
    nan_indices = np.isnan(values)

    # Create original grid points for interpolation
    grid_x, grid_z = np.meshgrid(x, z)
    points = np.array([grid_x.ravel(), grid_z.ravel()]).T
    values_flat = values.ravel()

    # Fill NaN values using nearest neighbor interpolation
    filled_values = griddata(
        points[~nan_indices.ravel()], 
        values_flat[~nan_indices.ravel()], 
        points, 
        method='nearest'
    ).reshape(values.shape)

    # Define new grid for interpolation
    new_depth = np.arange(0, max_z, delta_z)
    new_distance = np.arange(0, max_x, delta_x)
    new_grid_distance, new_grid_depth = np.meshgrid(new_distance, new_depth)
    new_points = np.array([new_grid_distance.ravel(), new_grid_depth.ravel()]).T

    # Perform linear interpolation on the new grid
    new_values = griddata(points, filled_values.ravel(), new_points, method='linear')

    # Reshape the interpolated values to match the new grid
    new_values_reshaped = new_values.reshape(new_grid_distance.shape)

    return new_values_reshaped, new_distance, new_depth


def filtrage_z(dz, Lcz, f):
    """
    Applies a low-pass Butterworth filter along the vertical (z) dimension of the dataset 
    to smooth the data and remove high-frequency noise.

    Parameters:
    ----------
    dz : float
        Grid spacing in the z-direction.
    Lcz : float
        Cutoff wavelength for filtering (higher values mean more smoothing).
    f : np.ndarray
        2D array of data to be filtered.

    Returns:
    -------
    f_smooth : np.ndarray
        Smoothed 2D array after applying the vertical filter.
    """

    f_smooth = np.full(f.shape, np.nan)  # Initialize with NaNs

    # Compute the sampling frequency and Nyquist frequency
    f_s = 1 / abs(dz)
    f_nyq = f_s / 2
    f_c = 1 / Lcz  # Cutoff frequency
    order = 4  # Filter order

    # Design a low-pass Butterworth filter
    b, a = signal.butter(order, f_c / f_nyq)

    for i in range(f.shape[0]):  # Loop over depth levels
        ind = np.where(~np.isnan(f[i, :]))[0]  # Find valid (non-NaN) indices
        if len(ind) > 15:
            # Remove the mean, filter, then restore the mean
            f_smooth[i, ind] = signal.filtfilt(b, a, signal.detrend(f[i, ind]))
            f_smooth[i, ind] += f[i, ind] - signal.detrend(f[i, ind])

    return f_smooth


def filtrage_x(dx, Lcx, f):
    """
    Applies a low-pass Butterworth filter along the horizontal (x) dimension of the dataset 
    to smooth the data and remove high-frequency noise.

    Parameters:
    ----------
    dx : float
        Grid spacing in the x-direction.
    Lcx : float
        Cutoff wavelength for filtering (higher values mean more smoothing).
    f : np.ndarray
        2D array of data to be filtered.

    Returns:
    -------
    f_smooth : np.ndarray
        Smoothed 2D array after applying the horizontal filter.
    """

    f_smooth = np.full(f.shape, np.nan)  # Initialize with NaNs

    # Compute the sampling frequency and Nyquist frequency
    f_s = 1 / abs(dx)
    f_nyq = f_s / 2
    f_c = 1 / Lcx  # Cutoff frequency
    order = 4  # Filter order

    # Design a low-pass Butterworth filter
    b, a = signal.butter(order, f_c / f_nyq)

    for j in range(f.shape[1]):  # Loop over horizontal positions
        ind = np.where(~np.isnan(f[:, j]))[0]  # Find valid (non-NaN) indices
        if len(ind) > 15:
            # Remove the mean, filter, then restore the mean
            f_smooth[ind, j] = signal.filtfilt(b, a, signal.detrend(f[ind, j]))
            f_smooth[ind, j] += f[ind, j] - signal.detrend(f[ind, j])

    return f_smooth