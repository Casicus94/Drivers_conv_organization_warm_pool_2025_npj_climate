import xarray as xr
import pandas as pd
import numpy as np
from xeofs.xarray import EOF, Rotator
import pdb


def is_sea(month,ini,fin):
    return (month >= ini) & (month <= fin)

user = '/home/netapp-clima/users/acasallas/ERA5_real/Area_2-9_135-145/'
scra = '/home/netapp-clima/scratch/acasallas/'

ds_sst = xr.open_dataset(user+'SST_area.nc')
ds_tcwv = xr.open_dataset(user+'TCWV_area.nc')
ds_u = xr.open_dataset(user+'Zonal_850_area.nc')
ds_v = xr.open_dataset(user+'Meridional_850_area.nc')
ds_u = ds_u.u.loc['2016-01-01 00:00:00':'2019-12-31 23:00:00']
ds_v = ds_v.v.loc['2016-01-01 00:00:00':'2019-12-31 23:00:00']

ds_sst = ds_sst.sst.groupby('time.month') - ds_sst.sst.groupby('time.month').mean()
ds_tcwv = ds_tcwv.tcwv.groupby('time.month') - ds_tcwv.tcwv.groupby('time.month').mean()

season = ['all','DJF','MAM','JJA','SON']

ini = [0,12,3,6,9]
fin = [0,2,5,8,11]
n_pc_es = 5

for i,sea in enumerate(season):
    print('We are calculating: '+sea)
    if sea == 'all':
        sst = ds_sst
        tcwv = ds_tcwv
        uu = ds_u
        vv = ds_v
    elif sea == 'DJF': 
        tmp_sst = []
        tmp_tcwv = []
        tmp_u = []
        tmp_v = []
        tmp_sst.append(ds_sst.sel(time=is_sea(ds_sst['time.month'],0,2)))
        tmp_sst.append(ds_sst.sel(time=is_sea(ds_sst['time.month'],12,12))) 
        tmp_tcwv.append(ds_tcwv.sel(time=is_sea(ds_tcwv['time.month'],0,2)))
        tmp_tcwv.append(ds_tcwv.sel(time=is_sea(ds_tcwv['time.month'],12,12)))
        tmp_u.append(ds_u.sel(time=is_sea(ds_u['time.month'],0,2)))
        tmp_u.append(ds_u.sel(time=is_sea(ds_u['time.month'],12,12)))
        tmp_v.append(ds_v.sel(time=is_sea(ds_v['time.month'],0,2)))
        tmp_v.append(ds_v.sel(time=is_sea(ds_v['time.month'],12,12)))
        con_sst = xr.concat(tmp_sst,dim='time')
        con_tcwv = xr.concat(tmp_tcwv,dim='time')
        con_u = xr.concat(tmp_u,dim='time')
        con_v = xr.concat(tmp_v,dim='time') 
        sst = con_sst
        tcwv = con_tcwv
        uu = con_u
        vv = con_v
    else:
        tmp_sst = ds_sst.sel(time=is_sea(ds_sst['time.month'],ini[i],fin[i]))
        sst = tmp_sst
        tmp_tcwv = ds_tcwv.sel(time=is_sea(ds_tcwv['time.month'],ini[i],fin[i]))
        tcwv = tmp_tcwv 
        tmp_u = ds_u.sel(time=is_sea(ds_u['time.month'],ini[i],fin[i]))
        uu = tmp_u
        tmp_v = ds_v.sel(time=is_sea(ds_v['time.month'],ini[i],fin[i]))
        vv = tmp_v
    ### Calculation
    #mpca = EOF([sst,tcwv,uu[:,0,:,:],vv[:,0,:,:]], dim=['time'], norm=False, weights='coslat')
    mpca = EOF([sst,tcwv], dim=['time'], norm=False, weights='coslat')
    mpca.solve()
    rot = Rotator(mpca, n_rot=50)
    reofs = rot.eofs()
    rpcs = rot.pcs()
    ##
    data_array0 = xr.DataArray(reofs[0], coords={'lat':ds_tcwv.latitude.values,'lon':ds_tcwv.longitude.values,'eofs':np.arange(0,50)}, dims=['lat','lon','eofs'])
    data_array1 = xr.DataArray(reofs[1], coords={'lat':ds_tcwv.latitude.values,'lon':ds_tcwv.longitude.values,'eofs':np.arange(0,50)}, dims=['lat','lon','eofs'])
    #data_array2 = xr.DataArray(reofs[2], coords={'lat':ds_tcwv.latitude.values,'lon':ds_tcwv.longitude.values,'eofs':np.arange(0,50)}, dims=['lat','lon','eofs'])
    #data_array3 = xr.DataArray(reofs[3], coords={'lat':ds_tcwv.latitude.values,'lon':ds_tcwv.longitude.values,'eofs':np.arange(0,50)}, dims=['lat','lon','eofs'])

    data_array0.to_netcdf(scra+'REOFs_SST_'+sea+'.nc')
    data_array1.to_netcdf(scra+'REOFs_TCWV_'+sea+'.nc')
    #data_array2.to_netcdf(scra+'REOFs_Zonal_'+sea+'.nc')
    #data_array3.to_netcdf(scra+'REOFs_Meridional_'+sea+'.nc')
    ###
    df_pcs = pd.DataFrame(data=rpcs)
    df_pcs.to_csv(scra+'PCAs_Rot_'+sea+'.nc') 
    ###
