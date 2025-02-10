import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset,num2date,date2num
from os.path import exists
import os
import pdb

path = '/home/esp-shared-a/Observations/SST/GHRSST/WPACIFIC/v200_nc4_normal_std/'
scra = '/home/netapp-clima/scratch/acasallas/'

years = np.arange(2016,2020)

for year in years:
    print('-----> Starting year: '+str(year))
    ini_dates = [str(year)+'-01-01 00:00:00', str(year)+'-04-01 00:00:00', str(year)+'-07-01 00:00:00', str(year)+'-10-01 00:00:00']
    end_date = [str(year)+'-03-31 23:00:00', str(year)+'-06-30 23:00:00', str(year)+'-09-30 23:00:00', str(year)+'-12-31 23:00:00']
    
    for j,ini_date in enumerate(ini_dates):
        horas = [d.strftime('%Y%m%d%H%M%S') for d in pd.date_range(ini_date,end_date[j], freq='1H')]
        datasets = np.ones((len(horas),851,2001))*np.nan
        count = 0
        for i,hour in enumerate(horas):
            print(hour)
            if exists(path+hour[0:6]+'/'+hour[6:8]+'/'+hour+'-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v2.0-v02.0-fv01.0.nc'): 
                ds = xr.open_dataset(path+hour[0:6]+'/'+hour[6:8]+'/'+hour+'-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v2.0-v02.0-fv01.0.nc')
                ds = ds['sea_surface_temperature'].sel(lat=slice(12,-5), lon = slice(130,170))
                datasets[i,:,:] = np.array(ds[0,:,:])
            else:
                print('No file for: '+hour)
        print('Creating file')
        unout = 'days since '+ini_date
        dates = [d for d in pd.date_range(ini_date,end_date[j], freq='1H')]
        datesout = date2num(dates,unout)
        ncout = Dataset(scra+'SST_Him_'+ini_date[0:7]+'.nc','w','NETCDF4')
        ncout.createDimension('longitude',2001)
        ncout.createDimension('latitude',851)
        ncout.createDimension('time',len(dates))
        lonvar=ncout.createVariable('longitude','float32',('longitude'));lonvar.setncattr('units','degrees_east');lonvar[:]=np.linspace(130,170,2001);
        latvar=ncout.createVariable('latitude','float32',('latitude'));latvar.setncattr('units','degrees_north');latvar[:]=np.array(np.linspace(12,-5,851));
        timevar = ncout.createVariable('time','float32',('time'));timevar.setncattr('units',unout);timevar[:]=datesout;

        myvar = ncout.createVariable('sea_surface_temperature','float32',('time','latitude','longitude'));myvar[:] = datasets;
        ncout.close();
        print('File created')
    print('CDO Merging year: '+str(year))
    os.system('cdo -b F32 mergetime '+scra+'SST_Him_'+str(year)+'* '+scra+'SST_Him_'+str(year)+'.nc')
    os.system('rm '+scra+'SST_Him_'+str(year)+'-*')
os.system('cdo -O -b F32 mergetime '+scra+'SST_Him_201* '+scra+'SST_Him_all.nc')

print('Complete')
