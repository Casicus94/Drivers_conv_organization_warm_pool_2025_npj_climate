import xarray as xr
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy import stats
import Casicus as casi
import pandas as pd

#Directories
user = '/home/tompkins-archive2/acasallas/ERA5_real/'
mdir = '/home/tompkins-archive/acasallas/Real1_run/Observations/'
data = '/home/tompkins-archive/acasallas/Data_to_plot/'

### Dates
init_year = '2016'
init_mon = '01'
init_day = '01'

end_year   = '2019'
end_mon   = '12'
end_day   = '31'

days = [d.strftime('%Y-%m-%d %H:%M:%S') for d in pd.date_range(init_year+init_mon+init_day,end_year+end_mon+end_day, freq='1H')]

###Area
latitud = '2-9'
longitud = '245-255'

###First we need to preprocess
print('Starting Hovmuller script!')
time1 = 0*24
time2 = 1460*24
ntime=time2-time1+1 # 240 quick hack should be read from the netcdf file really...
qv_press=1100. # low press qv integration hPA

dtime=3600 #This depends on the ouput frequency
difftime=str(dtime)

fstr=difftime+'-'+str(time1)+'-'+str(time2)

#Domain and splits
nx= 57 #29 #57 #15
ny= 81 #41 #81 #21

npts=nx*ny
hovchunk=81  #81 if 81 lon   # averaging length for chunking, power of 2, larger=smoother.
nhov=int(npts/hovchunk)
hovx=100*(np.arange(nhov)+0.5)/nhov

#Start matrixes 
times=np.array(range(ntime)) # start from day 1 for spin up
#slices=times+time1
hovtimes=(times+time1)*dtime/86400.
sst_hov=np.zeros([1,ntime,nhov])
tcwv_hov=np.zeros([1,ntime,nhov])

print('Area_'+latitud+'_'+longitud)
sst_ds = xr.open_dataset(user+'Area_'+latitud+'_'+longitud+'/SST_area.nc')
tcwv_ds = xr.open_dataset(user+'Area_'+latitud+'_'+longitud+'/TCWV_area.nc')
sst_ds = sst_ds.sel(time=slice('2016-01-01 00:00:00','2019-12-31 23:00:00'))
tcwv_ds = tcwv_ds.sel(time=slice('2016-01-01 00:00:00','2019-12-31 23:00:00'))

#Calculate indexes
print('Calculating IQR')
IQR = np.array([])
variable = 'tcwv'
for i in range(1,ntime):
    CRH = tcwv_ds[variable][i,:,:] 
    IQR = np.append(IQR, np.quantile(CRH,0.75) - np.quantile(CRH,0.25))

for itime in times:
    itime1 = itime+time1
    slice=itime1+time1
    print("slice ",itime1,slice)
    tcwv = np.array(tcwv_ds.variables['tcwv'][slice,:,:])
    sst=  np.array(sst_ds.variables["sst"][slice,:,:])
    print("min max tcwv ",np.min(tcwv),np.max(tcwv))
    # sort the tcwv
    isort=np.argsort(tcwv.flatten()) # sort tcwv
    sstpert=sst.flatten()-np.mean(sst)
    sst_hov[0,itime,:]=np.mean(sstpert[isort].reshape(-1,hovchunk),axis=1)
    tcwv_hov[0,itime,:]=np.mean(tcwv.flatten()[isort].reshape(-1,hovchunk),axis=1)
#   
# calculate slopes
#
x = np.array(hovx)
slopes = []
for i in range(1, ntime):
    y = np.array(sst_hov[0][i])
    xm = np.ma.masked_array(x,mask=np.isnan(y)).compressed()
    ym = np.ma.masked_array(y,mask=np.isnan(y)).compressed()
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(xm, ym)
        slopes.append(slope)  
    except:
        print("Error: ",i)
        slopes.append(np.nan)
# negative standard deviation 
slopes_smth=pd.Series(slopes).rolling(window=24,min_periods=1,center=True).mean()
slo_data = pd.DataFrame(np.column_stack((days[:-1],np.array(slopes_smth),np.array(IQR))), columns=['Dates','Slope','IQR'])
slo_data.to_csv(data+'Slopes_data_'+latitud+'_'+longitud+'_ERA5.csv',index=False)

