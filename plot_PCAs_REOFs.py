import pandas as pd
import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *

user = '/home/tompkins-archive2/acasallas/ERA5_real/Area_2-9_135-145/'
scra = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive2/acasallas/REOFs/'

df_slope = pd.read_csv(user+'Slopes_data_2-9_135-145_ERA5.csv')
df_reof = pd.read_csv(mdir+'PCAs_Rot_all.nc')
df_slope['PCS'] = (df_reof['2'][0:35040]+df_reof['1'][0:35040]+df_reof['0'][0:35040])*(-1)

df_slope.set_index('Dates', inplace=True)
df_slope1 = df_slope.loc['2017-03-15 00:00:00':'2017-05-15 23:00:00']

dates_MAM = [d.strftime('%Y-%m-%d') for d in pd.date_range('2017-03-15 00:00:00','2017-05-15 23:00:00', freq='6d')]
tval_MAM = pd.date_range('2017-03-15 00:00:00','2017-05-15 23:00:00', freq='1H')

fig = plt.figure(figsize=(12,4))
gs = GridSpec(1,1,left = 0.1, right = 0.96, hspace=0.15, wspace=0.15, top = 0.95, bottom = 0.08)
ax = plt.subplot(gs[0])
plt.plot(tval_MAM,df_slope1['Slope'], color = 'blue', label = 'Slope', alpha = 0.8)
plt.plot(tval_MAM,df_slope1['PCS'], color = 'purple', label = 'Sum of PCs 1-2-3', alpha = 0.8)
plt.axhline(df_slope['Slope'].mean()+0.5*df_slope['Slope'].std(), color = 'blue', linestyle = '--')
plt.axhline(df_slope['Slope'].mean()-0.5*df_slope['Slope'].std(), color = 'blue', linestyle = '--')
plt.axhline(df_slope['PCS'].mean()+0.5*df_slope['PCS'].std(), color = 'purple', linestyle = '--')
plt.axhline(df_slope['PCS'].mean()-0.5*df_slope['PCS'].std(), color = 'purple', linestyle = '--')
plt.legend(frameon=False)
plt.text(tval_MAM[0],-0.02,'$R_{Pearson}$=0.87 for all boreal winter/spring months', fontsize = 14)
plt.text(tval_MAM[0],-0.015,'$R_{Pearson}$=0.72 for the entire period', fontsize = 14)
ax.set_xticks(dates_MAM)
ax.set_xticklabels(dates_MAM)
plt.savefig(scra+'PCA1-3_Slope.jpg',dpi=300, bbox_inches = 'tight')
plt.savefig(scra+'PCA1-3_Slope.pdf',dpi=300, bbox_inches = 'tight')
plt.show()
plt.close()

pdb.set_trace()
##############################
df_pcs = (df_reof['0']+df_reof['1']+df_reof['2'])*(-1)
dates = df_slope.Dates

df = pd.DataFrame(np.column_stack((dates,df_pcs[0:len(dates)])))
df.to_csv(user+'REOFs_slope.csv', index=False)


