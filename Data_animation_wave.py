import xarray as xr
import numpy as np
import pdb
import datetime as dt
import os
import pandas as pd

def sum_them(ds_tcwv,ds_olr,ds_mu,ds_mv,ds_sst,dates):
    date_olr = dates
    data_olr = []
    for i,date in enumerate(date_olr):
         try:
             data_olr.append(ds_olr.ttr.loc[date])
         except:
             print('Date not available')
    weights = np.cos(np.deg2rad(ds_olr['latitude']))
    tmp_olr = xr.concat(data_olr, dim = 'time')
    olrw = tmp_olr.weighted(weights)
    olrtw = ds_olr.groupby('time.season').mean()
    olr_rev = olrw.mean(dim='time') - olrtw.ttr[2]
    data_tcwv = []
    data_mu = []
    data_mv = []
    data_sst = []
    for i,date in enumerate(dates):
        try:
            data_tcwv.append(ds_tcwv.tcwv.loc[date])
            data_mu.append(ds_mu['p71.162'].loc[date])
            data_mv.append(ds_mv['p72.162'].loc[date])
            data_sst.append(ds_sst.sst.loc[date])
        except:
            print('Date not available')
    ###Means
    #TCWV
    tcwv_tmp = xr.concat(data_tcwv, dim='time')
    tcwvmean = ds_tcwv.groupby('time.season').mean()
    tcwv_rev = tcwv_tmp.mean(dim='time')-tcwvmean.tcwv[2]
    #SST
    sst_tmp = xr.concat(data_sst, dim='time')
    sst_rev = sst_tmp.mean(dim='time')
    #Moisture flux eastwar and northward
    mu_tmp = xr.concat(data_mu, dim='time')
    mu_rev = mu_tmp.mean(dim='time')
    mv_tmp = xr.concat(data_mv, dim='time')
    mv_rev = mv_tmp.mean(dim='time')
    return(tcwv_rev,olr_rev,mu_rev,mv_rev,sst_rev)

def sel_season(data,sta,med,end):
    tmp = data[data.index.month.isin([sta,med,end])]
    new_df = tmp.index.strftime('%Y-%m-%d %H:%M:%S').values
    return(new_df)

def means_cal(ds_tcwv,ds_olr,ds_mu,ds_mv,ds_sst,latitude,longitude,mdir,season,sta,med,end):
    scra = '/home/netapp-clima/scratch/acasallas/Composites_'+latitude+'/'
    tipo = '' #'' for normal slope, or REOFs for PCAs
    dates_rev = pd.read_csv(mdir+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/Dates_'+latitude+'_'+longitude+'_reversals'+tipo+'.csv', parse_dates = True, index_col = 0)
    dates_org = pd.read_csv(mdir+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/Dates_'+latitude+'_'+longitude+'_organised'+tipo+'.csv', parse_dates = True, index_col = 0)
    
    if season == 'all':
        dat_rev = dates_rev.index.strftime('%Y-%m-%d %H:%M:%S').values
        dat_org = dates_org.index.strftime('%Y-%m-%d %H:%M:%S').values
    elif season == 'DJF' or season == 'MAM' or season == 'JJA' or season == 'SON':
        dat_rev = sel_season(dates_rev,sta,med,end) 
        dat_org = sel_season(dates_org,sta,med,end)
    else:
        tmp_rev = dates_rev[dates_rev.index.month.isin([3,4,5])]
        tmp_org = dates_org[dates_org.index.month.isin([3,4,5])]
        dat_rev = (tmp_rev.index -  pd.Timedelta(days=float(season))).strftime('%Y-%m-%d %H:%M:%S').values
        dat_org = (tmp_org.index - pd.Timedelta(days=float(season))).strftime('%Y-%m-%d %H:%M:%S').values
    
    tcwv_rev,olr_rev,mu_rev,mv_rev,sst_rev = sum_them(ds_tcwv,ds_olr,ds_mu,ds_mv,ds_sst,dat_rev)
    tcwv_org,olr_org,mu_org,mv_org,sst_org = sum_them(ds_tcwv,ds_olr,ds_mu,ds_mv,ds_sst,dat_org)
    nam_rev = ['TCWV_rev','OLR_rev','mu_rev','mv_rev','SST_rev']
    nam_org = ['TCWV_org','OLR_org','mu_org','mv_org','SST_org'] 
    revers = [tcwv_rev,olr_rev,mu_rev,mv_rev,sst_rev]
    org = [tcwv_org,olr_org,mu_org,mv_org,sst_org]
    for i,rev in enumerate(revers):
        rev.to_netcdf(scra+nam_rev[i]+'_'+season+tipo+'_'+latitude+'.nc') 
        org[i].to_netcdf(scra+nam_org[i]+'_'+season+tipo+'_'+latitude+'.nc')

    return(tcwv_rev,olr_rev,mu_rev,mv_rev,sst_rev, tcwv_org,olr_org,mu_org,mv_org,sst_org)

mdir = '/home/tompkins-archive/acasallas/'
user = '/home/tompkins-archive2/acasallas/ERA5_real/'

print('Reading TCWV')
ds_tcwv = xr.open_dataset(mdir+'ERA5/TCWV_all.nc')
print('Reading OLR')
ds_olr = xr.open_dataset(user+'TOA_net_LW_all.nc')
print('Reading SST')
ds_sst = xr.open_dataset(mdir+'ERA5/SST_all.nc')
print('Reading Moisture flux')
ds_mu = xr.open_dataset(user+'UWVF_integral_all.nc')
ds_mv = xr.open_dataset(user+'VWVF_integral_all.nc')

latitude = '3-10'
longitude = '147-157'
scra = '/home/netapp-clima/scratch/acasallas/Composites_'+latitude+'/'

print('Running '+latitude)
season = ['all','DJF','MAM','JJA','SON']
sta = [0,12,3,6,9]
med = [0,1,4,7,10]
end = [0,2,5,8,11]

sea_ind = False
if sea_ind == True:
    for i,sea in enumerate(season):
        print('We are calculating '+sea)
        tcwv_rev,olr_rev,sst_rev, tcwv_org,olr_org,sst_org = means_cal(ds_tcwv,ds_olr,ds_sst,latitude,longitude,mdir,sea,sta[i],med[i],end[i])

seas = []
for i in np.arange(-9,9.00001,1/24):
    seas.append(str(i))

#seas = ['9','8.5','8','7.5','7','6','5','4.5','4','3.5','3','2','1','0','-1','-2','-3']
lag = True
if lag == True:
    print('Calculating lag')
    for sea in seas:
        print('Calculating time: '+sea)
        if not os.path.exists(scra+'TCWV_rev_'+sea+'_'+latitude+'.nc'): 
            tcwv_rev,olr_rev,mu_rev,mv_rev,sst_rev, tcwv_org,olr_org,mu_org,mv_org,sst_org = means_cal(ds_tcwv,ds_olr,ds_mu,ds_mv,ds_sst,latitude,longitude,mdir,sea,0,0,0)


