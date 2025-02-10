import numpy as np
import h5py
import xarray as xr
from netCDF4 import Dataset,num2date,date2num
import pandas as pd 
import pdb
import os 

gpm_path = '/home/esp-shared-a/Observations/GPM/IMERG/FINAL/V06B/'
scra = '/home/netapp-clima/scratch/acasallas/'
pre_file = '3B-HHR.MS.MRG.3IMERG.'

horas = []
for i in np.arange(0,1411,30):
    if i == 0:
        horas.append('00'+str(i)+'0')
    elif i < 100:
        horas.append('00'+str(i))
    elif i < 1000:
        horas.append('0'+str(i))
    else:
        horas.append(str(i))

#There must be better ways to do this!
index = ['002','005','012','015','022','025','032','035','042','045','052','055','062',
         '065','072','075','082','085','092','095','102','105','112','115','122','125',
         '132','135','142','145','152','155','162','165','172','175','182','185','192',
         '195','202','205','212','215','222','225','232','235']

day_hrs = [d.strftime('%H%M%S') for d in pd.date_range('2015-01-01 00:00:00','2015-01-01 23:30:00',freq='0.5H')]
years = np.arange(2015,2021)
tot_mat = np.ones((len(years),1800,3600))*np.nan

#I have to do it year by year due to the ram
for iyear,year in enumerate(years):
    print('-----> Starting year: '+str(year))
    sta_date = str(year)+'-01-01 00:00:00'
    end_date = str(year)+'-12-31 23:30:00'
    dat_for = [d.strftime('%Y%m%d') for d in pd.date_range(sta_date,end_date,freq='1D')]
    year_mat = np.ones((len(dat_for),1800,3600))*np.nan
    for idate,fdate in enumerate(dat_for):
        print('------- '+fdate+' -------')
        day_mat = np.ones((len(day_hrs[:-1]),1800,3600))*np.nan
        for i,hrs in enumerate(day_hrs[:-1]):
            print('##### '+hrs+' #####')
            try: 
                f = h5py.File(gpm_path+fdate[0:4]+'/'+fdate[4:6]+'/'+fdate[6:8]+'/'+pre_file+fdate+'-S'+hrs+'-E'+index[i]+'959.'+horas[i]+'.V06B.HDF5')
                precip = f['Grid/IRprecipitation'][0][:][:]
                del(f)
                precip = np.transpose(precip)
                day_mat[i,:,:] = precip
                del(precip)
            except:
                print('File Missing')
        day_mat[day_mat<0] = np.nan
        sum_day = np.nansum(day_mat, axis = 0)
        sum_day[sum_day==0] = np.nan
        del(day_mat)
        year_mat[idate,:,:] = sum_day
        del(sum_day)
    #tot_mat[iyear,:,:] = np.nanmean(year_mat, axis = 0)
    #del(year_mat)
    year_mat[year_mat==np.nan] = 0
    print('Creating file')
    unout = 'days since '+sta_date
    datess = [d for d in pd.date_range(sta_date,end_date, freq='1D')]
    datesout = date2num(datess,unout)
    ncout = Dataset(scra+'Precip_GPM_daily_'+str(year)+'.nc','w','NETCDF4')
    ncout.createDimension('longitude',3600)
    ncout.createDimension('latitude',1800)
    ncout.createDimension('time',len(dat_for))
    lonvar=ncout.createVariable('longitude','float32',('longitude'));lonvar.setncattr('units','degrees_east');lonvar[:]=np.arange(-179.95,179.951,0.1);
    latvar=ncout.createVariable('latitude','float32',('latitude'));latvar.setncattr('units','degrees_north');latvar[:]=np.arange(-89.95,89.951,0.1);
    timevar = ncout.createVariable('time','float32',('time'));timevar.setncattr('units',unout);timevar[:]=datesout;
    myvar = ncout.createVariable('Precip','float32',('time','latitude','longitude'));myvar[:] = year_mat;
    ncout.close();
    print('File created')
print('Merging files')
os.system('cdo -b F32 mergetime '+scra+'Precip_GPM_daily_20* '+scra+'Precip_GPM_daily_all.nc')
print('Complete')
#### ds.filna(0)
