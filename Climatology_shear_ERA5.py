mport xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.gridspec as gridspec
from pylab import *
import Casicus as casi
import pandas as pd

#Directories
share = '/home/esp-shared-a/Observations/CERES/SYN1deg/'
scr = '/home/netapp-clima/scratch/acasallas/wrf/Plots/'
rea = '/home/netapp-clima-shared/MODELS/REANALYSIS/ERA5/'
mdir = '/home/tompkins-archive/acasallas/Real1_run/Observations/'

#Read data
dsu = xr.open_dataset(mdir+'u_component_of_wind_15-20.nc')
dsv = xr.open_dataset(mdir+'v_component_of_wind_15-20.nc')

#Shear 1 = 700 hPa and 2 = 850 hPa
du = (dsu.u[:,1,:,:] - dsu.u[:,2,:,:])
dv = (dsv.v[:,1,:,:] - dsv.v[:,2,:,:])

#Monthly mean
du_mon = du.groupby("time.month").mean()
dv_mon = dv.groupby("time.month").mean()
#Daily Mean
du_yday = du.groupby("time.dayofyear").mean()
dv_yday = dv.groupby("time.dayofyear").mean()
du_day = du.resample(time='d').mean()
dv_day = dv.resample(time='d').mean()

#fldmean
du_monmean = du_mon.mean(dim=['latitude','longitude'])
dv_monmean = dv_mon.mean(dim=['latitude','longitude'])
du_ydaymean = du_yday.mean(dim=['latitude','longitude'])
dv_ydaymean = dv_yday.mean(dim=['latitude','longitude'])
du_daymean = du_day.mean(dim=['latitude','longitude'])
dv_daymean = dv_day.mean(dim=['latitude','longitude'])

years = np.arange(2016,2021)
plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(2,1, left=0.1, right=0.975, hspace=0.25, wspace=0.0, top=0.9, bottom=0.1)

### Zonal shear
ax = subplot(gs[0])

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
for i,year in enumerate(years):
    day = du_daymean.sel(time=slice(str(year)+'-01-01',str(year)+'-12-31'))
    plt.plot(np.linspace(1,13,len(day)),day)
plt.plot(np.linspace(1,13,len(du_ydaymean)),du_ydaymean, color = 'k', linewidth = 1.3)
plt.axhline(0, linestyle=':', color = 'k')
plt.xticks(np.arange(1,13))
plt.ylabel('Zonal shear (m s$^{-1}$ km$^{-1}$)')
plt.xlim(1,13)
plt.xlabel('Time (Months)')

### Meridional Shear
ax = subplot(gs[1])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
for i,year in enumerate(years):
    day = dv_daymean.sel(time=slice(str(year)+'-01-01',str(year)+'-12-31'))
    plt.plot(np.linspace(1,13,len(day)),day)
plt.plot(np.linspace(1,13,len(dv_ydaymean)),dv_ydaymean, color = 'k', linewidth = 1.3)
plt.axhline(0, linestyle=':', color = 'k')
plt.xticks(np.arange(1,13))
plt.ylabel('Meriodional shear (m s$^{-1}$ km$^{-1}$)')
plt.xlim(1,13)
plt.xlabel('Time (Months)')

plt.savefig(scr+'Climatology_shear.jpg',  bbox_inches = 'tight')
plt.show()

