import pandas as pd
import numpy as np
import xarray as xr
from sklearn.metrics import mean_squared_error
import pdb 

mdir = '/home/tompkins-archive/acasallas/Data_to_plot/'

ds = xr.open_dataset(mdir+'slopes_adrian.nc')
df = pd.read_csv(mdir+'Slopes_data_2-9_135-145_ERA5.csv', parse_dates = True, index_col = 0)
df_cut = df.loc['2017-03-14 15:00:00':'2017-05-09 14:00:00']

rms = mean_squared_error(ds.slope.values, df_cut['Slope'], squared=True)
print('The RMSE is: '+str(rms))
print('The correlation is: '+str(np.corrcoef(ds.slope.values, df_cut['Slope'])[0,1]))

pdb.set_trace()
