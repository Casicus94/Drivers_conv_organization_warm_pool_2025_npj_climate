import pandas as pd
import xarray as xr
import numpy as np
import pdb 

### Dates
ini_date = '2016-10-01 00:00:00'
end_date = '2020-01-01 00:00:00'

horas = [d.strftime('%Y-%m-%d %H:%M:%S') for d in pd.date_range(ini_date,end_date, freq='1D')]

mdir = '/home/tompkins-archive2/acasallas/Observations/'
###Area
area = '2' #2-9 = '' 3-10 = '1' 3s-4n = '2'
tcwv_ds = xr.open_dataset(mdir+'Mimic_area'+area+'.nc')

#Calculate indexes
print('Calculating IQR')
IQR = np.array([])
variable = 'tpwGrid'
ntime = len(horas)
for i in range(0,ntime):
    CRH = tcwv_ds[variable][i,:,:]
    IQR = np.append(IQR, np.quantile(CRH,0.75) - np.quantile(CRH,0.25))

data = pd.DataFrame(np.column_stack((horas,np.array(IQR))), columns=['Dates','IQR_pw'])
if area == '':
    data.to_csv(mdir+'TPW_IQR_area_2-9_135-145.csv', index = False)
elif area == '1':
    data.to_csv(mdir+'TPW_IQR_area_3-10_147-157.csv', index = False)
elif area == '2':
    data.to_csv(mdir+'TPW_IQR_area_3S-4N_156-166.csv', index = False)

