import xarray as xr
import numpy as np
import pdb
import datetime as dt
import os
from os.path import exists
import pandas as pd

def sum_them(ds_tcwv,ds_olr,ds_u,ds_v,ds_div,ds_md,ds_mu,ds_mv,ds_sst,ds_ross,dates):
    date_ross = []
    for date in dates:
        date_ross.append(date[0:10])
    date_ross = np.unique(date_ross)
    data_ross = []   
    for i,date in enumerate(date_ross):
         try:
             data_ross.append(ds_ross.Source_term.loc[date])
         except:
             print('Date not available')
    #Source Term
    sot_tmp = xr.concat(data_ross, dim='time')
    sot_rev = sot_tmp.mean(dim='time')

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
    data_u = []
    data_v = []
    d_d850 = []
    d_d250 = []
    data_md = []
    data_mu = []
    data_mv = []
    data_sst = []
    for i,date in enumerate(dates):
        try:
            data_tcwv.append(ds_tcwv.tcwv.loc[date])
            data_u.append(ds_u.u.loc[date])
            data_v.append(ds_v.v.loc[date])
            d_d850.append(ds_div.d.loc[date][1])
            d_d250.append(ds_div.d.loc[date][0])
            data_md.append(ds_md['p84.162'].loc[date])
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
    #Wind
    u_tmp = xr.concat(data_u, dim='time')
    u_rev = u_tmp.mean(dim='time')
    v_tmp = xr.concat(data_v, dim='time')
    v_rev = v_tmp.mean(dim='time')
    #Divergence
    div850_tmp = xr.concat(d_d850, dim='time')
    divmean = ds_div.groupby('time.season').mean()
    #d850_rev = div850_tmp.mean(dim='time')-divmean.d[2][1]
    d850_rev = div850_tmp.mean(dim='time')
    div250_tmp = xr.concat(d_d250, dim='time')
    #d250_rev = div250_tmp.mean(dim='time')-divmean.d[2][0]
    d250_rev = div250_tmp.mean(dim='time')
    #Moisture divergence flux
    md_tmp = xr.concat(data_md, dim='time')
    mdmean = ds_md.groupby('time.season').mean()
    md_rev = md_tmp.mean(dim='time')-mdmean['p84.162'][2]
    #Moisture flux eastwar and northward
    mu_tmp = xr.concat(data_mu, dim='time')
    mu_rev = mu_tmp.mean(dim='time')
    mv_tmp = xr.concat(data_mv, dim='time')
    mv_rev = mv_tmp.mean(dim='time')
    return(tcwv_rev,olr_rev,u_rev,v_rev,d850_rev,d250_rev,md_rev,mu_rev,mv_rev,sst_rev,sot_rev)

def sel_season(data,sta,med,end):
    tmp = data[data.index.month.isin([sta,med,end])]
    new_df = tmp.index.strftime('%Y-%m-%d %H:%M:%S').values
    return(new_df)

def means_cal(ds_tcwv,ds_olr,ds_u,ds_v,ds_div,ds_md,ds_mu,ds_mv,ds_sst,ds_ross,latitude,longitude,mdir,season,sta,med,end):
    scra = '/home/netapp-clima/scratch/acasallas/Composites_3-10/'
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
    
    print('Calculating means/anomalies!')
    tcwv_rev,olr_rev,u_rev,v_rev,d850_rev,d250_rev,md_rev,mu_rev,mv_rev,sst_rev,sot_rev = sum_them(ds_tcwv,ds_olr,ds_u,ds_v,ds_div,ds_md,ds_mu,ds_mv,ds_sst,ds_ross,dat_rev)
    tcwv_org,olr_org,u_org,v_org,d850_org,d250_org,md_org,mu_org,mv_org,sst_org,sot_org = sum_them(ds_tcwv,ds_olr,ds_u,ds_v,ds_div,ds_md,ds_mu,ds_mv,ds_sst,ds_ross,dat_org)
    nam_rev = ['TCWV_rev','OLR_rev','u_rev','v_rev','d850_rev','d250_rev','md_rev','mu_rev','mv_rev','SST_rev','Ross_rev']
    nam_org = ['TCWV_org','OLR_org','u_org','v_org','d850_org','d250_org','md_org','mu_org','mv_org','SST_org','Ross_org'] 
    revers = [tcwv_rev,olr_rev,u_rev,v_rev,d850_rev,d250_rev,md_rev,mu_rev,mv_rev,sst_rev,sot_rev]
    org = [tcwv_org,olr_org,u_org,v_org,d850_org,d250_org,md_org,mu_org,mv_org,sst_org,sot_org]
    print('Creating files')
    for i,rev in enumerate(revers):
        if len(ds_ross.coords) == 3:
            rev = rev.drop_vars('season')
            org[i] = org[i].drop_vars('season')
        rev.to_netcdf(scra+nam_rev[i]+'_'+season+tipo+'_'+latitude+'.nc') 
        org[i].to_netcdf(scra+nam_org[i]+'_'+season+tipo+'_'+latitude+'.nc')

    return(tcwv_rev,olr_rev,u_rev,v_rev,d850_rev,d250_rev,md_rev,mu_rev,mv_rev,sst_rev,sot_rev, tcwv_org,olr_org,u_org,v_org,d850_org,d250_org,md_org,mu_org,mv_org,sst_org,sot_org)

scra = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/'
user = '/home/tompkins-archive2/acasallas/ERA5_real/'

print('Reading TCWV')
ds_tcwv = xr.open_dataset(mdir+'ERA5/TCWV_all.nc')
print('Reading OLR')
ds_olr = xr.open_dataset(user+'TOA_net_LW_all.nc')
print('Reading SST')
ds_sst = xr.open_dataset(mdir+'ERA5/SST_all.nc')
print('Reading wind')
ds_u = xr.open_dataset(mdir+'ERA5/Zonal_850.nc')
ds_v = xr.open_dataset(mdir+'ERA5/Meridional_850.nc')
print('Reading divergence')
ds_div = xr.open_dataset(user+'divergence_all.nc')
print('Reading Moist div')
ds_md = xr.open_dataset(user+'Vert_int_div_moist_flux_all.nc')
print('Reading Moisture flux')
ds_mu = xr.open_dataset(user+'UWVF_integral_all.nc')
ds_mv = xr.open_dataset(user+'VWVF_integral_all.nc')
print('Reading Rossby Source Term')
ds_ross = xr.open_dataset(user+'Rossby_wave.nc')

latitude = '3-10'
longitude = '147-157'

season = ['all','DJF','MAM','JJA','SON']
sta = [0,12,3,6,9]
med = [0,1,4,7,10]
end = [0,2,5,8,11]

sea_ind = False
if sea_ind == True:
    for i,sea in enumerate(season):
        print('We are calculating '+sea)
        tcwv_rev,olr_rev,u_rev,v_rev,sst_rev,sot_rev, tcwv_org,olr_org,u_org,v_org,sst_org,sot_org = means_cal(ds_tcwv,ds_olr,ds_u,ds_v,ds_sst,ds_ross,latitude,longitude,mdir,sea,sta[i],med[i],end[i])

#seas = []
#for i in np.arange(-9,9.00001,1/24):
#    seas.append(str(i))

seas = ['9','8.5','8','7.5','7','6','5','4.5','4','3.5','3','2','1','0','-1','-2','-3']

lag = True
if lag == True:
    print('Calculating lag')
    for sea in seas:
        tcwv_rev,olr_rev,u_rev,v_rev,d850_rev,d250_rev,md_rev,mu_rev,mv_rev,sst_rev,sot_rev, tcwv_org,olr_org,u_org,v_org,d850_org,d250_org,md_org,mu_org,mv_org,sst_org,sot_org = means_cal(ds_tcwv,ds_olr,ds_u,ds_v,ds_div,ds_md,ds_mu,ds_mv,ds_sst,ds_ross,latitude,longitude,mdir,sea,0,0,0)

