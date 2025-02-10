import xarray as xr
import pandas as pd
import numpy as np
from eofs.multivariate.standard import MultivariateEof
import pdb


def is_sea(month,ini,fin):
    return (month >= ini) & (month <= fin)

user = '/home/netapp-clima/users/acasallas/ERA5_real/Area_2-9_135-145/'
scra = '/home/netapp-clima/scratch/acasallas/'

ds_sst = xr.open_dataset(user+'SST_area.nc')
ds_tcwv = xr.open_dataset(user+'TCWV_area.nc')

ds_sst = ds_sst.sst.groupby('time.month') - ds_sst.sst.groupby('time.month').mean()
ds_tcwv = ds_tcwv.tcwv.groupby('time.month') - ds_tcwv.tcwv.groupby('time.month').mean()

season = ['all','DJF','MAM','JJA','SON']

ini = [0,12,3,6,9]
fin = [0,2,5,8,11]
n_pc_es = 5

for i,sea in enumerate(season):
    print('We are calculating: '+sea)
    if sea == 'all':
        sst = np.array(ds_sst)
        tcwv = np.array(ds_tcwv)
    elif sea == 'DJF': 
        tmp_sst = []
        tmp_tcwv = []
        tmp_sst.append(ds_sst.sel(time=is_sea(ds_sst['time.month'],0,2)))
        tmp_sst.append(ds_sst.sel(time=is_sea(ds_sst['time.month'],12,12))) 
        tmp_tcwv.append(ds_tcwv.sel(time=is_sea(ds_tcwv['time.month'],0,2)))
        tmp_tcwv.append(ds_tcwv.sel(time=is_sea(ds_tcwv['time.month'],12,12)))
        con_sst = xr.concat(tmp_sst,dim='time')
        con_tcwv = xr.concat(tmp_tcwv,dim='time')
        sst = np.array(con_sst)
        tcwv = np.array(con_tcwv)
    else:
        tmp_sst = ds_sst.sel(time=is_sea(ds_sst['time.month'],ini[i],fin[i]))
        sst = np.array(tmp_sst)
        tmp_tcwv = ds_tcwv.sel(time=is_sea(ds_tcwv['time.month'],ini[i],fin[i]))
        tcwv = np.array(tmp_tcwv)     
        
    ### Calculation
    solver = MultivariateEof([sst, tcwv])

    ###
    s_pcs = pd.DataFrame(data=solver.pcs(npcs=n_pc_es, pcscaling=0), columns = ['PC1','PC2','PC3','PC4','PC5'])
    s_evfs = pd.DataFrame(data=solver.varianceFraction(neigs=n_pc_es), columns = ['VarianceFraction'])

    s_pcs.to_csv(scra+'PCAs_SST_TCWV_'+sea+'.csv')
    s_evfs.to_csv(scra+'Var_frac_SST_TCWV_'+sea+'.csv')

    ###
    data_array_0 = xr.DataArray(solver.eofs(neofs=n_pc_es, eofscaling=0)[0], coords={'eofs':np.arange(0,n_pc_es),'lat':ds_tcwv.latitude.values,'lon':ds_tcwv.longitude.values}, dims=['eofs','lat','lon'])
    data_array_1 = xr.DataArray(solver.eofs(neofs=n_pc_es, eofscaling=0)[1], coords={'eofs':np.arange(0,n_pc_es),'lat':ds_tcwv.latitude.values,'lon':ds_tcwv.longitude.values}, dims=['eofs','lat','lon'])

    data_array_0.to_netcdf(scra+'EOFs_SST_'+sea+'.nc')
    data_array_1.to_netcdf(scra+'EOFs_TCWV_'+sea+'.nc')

