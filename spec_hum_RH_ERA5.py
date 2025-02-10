import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as GridSpec
from pylab import *
import xarray as xr
import pdb

scra = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'
latitude = '2-9'
longitude = '135-145'

variable = 'RH' # RH spec_humidity

print('Reading data')
data = xr.open_dataset(mdir+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/'+variable+'_area.nc')

### Maybe is better if we first divide by regions, so these are the region coordinates
stas = [9,8,7,6,5,4,3]
end = [8,7,6,5,4,3,2]

means = {}
seasons = {}
for i,sta in enumerate(stas):
    if variable=='RH':
        tmp = data.r.sel(latitude=slice(sta,end[i]))
    elif variable=='spec_humidity':
        tmp = data.q.sel(latitude=slice(sta,end[i]))
    sea_tmp = tmp.groupby('time.season').mean()
    means[str(sta)+'-'+str(end[i])] = (tmp.mean(dim=['longitude','latitude','time']).values)
    seasons[str(sta)+'-'+str(end[i])] = sea_tmp.mean(dim=['longitude','latitude'])
    del(tmp);del(sea_tmp)

z = [100,125,150,175,200,225,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000] 

if variable == 'RH':
    cte = 1
elif variable == 'spec_humidity':
    cte = 1000

plt.figure(figsize=(5,8))
for i,sta in enumerate(stas):
    plt.plot(means[str(sta)+'-'+str(end[i])]*cte,z, label=str(end[i])+'N-'+str(sta)+'N')
plt.ylim(1000,100)
plt.ylabel('Pressure')
if variable == 'RH':
    plt.xlabel('Relative Humidity (%)')
elif variable == 'spec_humidity':
    plt.xlabel('Specific Humidity (g/kg)')
plt.legend(frameon=False)
plt.savefig(scra+variable+'_vert.jpg', bbox_inches='tight', dpi=300)
plt.savefig(med+variable+'_vert.jpg', bbox_inches='tight', dpi=300)
plt.savefig(med+variable+'_vert.pdf', bbox_inches='tight', dpi=300)
#plt.show()
plt.close()

#### Seasons
temps = ['(a) DJF','(b) MAM','(c) JJA','(d) SON']
fig = plt.figure(figsize=(8,8.19))
gs = GridSpec(1,4,left = 0.09, right = 0.95, hspace=0.35, wspace=0.13, top = 0.88, bottom = 0.2)
posi = [0,2,1,3]
fs = 8

for j,sea in enumerate(temps):
    ax = plt.subplot(gs[j])
    for i,sta in enumerate(stas):
        plt.plot(seasons[str(sta)+'-'+str(end[i])][posi[j]]*cte,z, label=str(end[i])+'N-'+str(sta)+'N')
    plt.ylim(1000,100)
    if j == 0:
        plt.ylabel('Pressure (hPa)', fontsize = fs)
        plt.yticks(fontsize = fs)
    else:
        ax.yaxis.set_major_formatter(NullFormatter())
    if variable == 'RH':
        plt.xlabel('Relative Humidity (%)', fontsize = fs)
    elif variable == 'spec_humidity':
        plt.xlabel('Specific Humidity (g/kg)', fontsize = fs)
    if j == 1:
    	plt.legend(loc='upper center',frameon=False, ncol=len(stas), bbox_to_anchor=(1.05,1.08), fontsize = fs)
    plt.title(sea,fontweight='bold', loc = 'left', fontsize = fs)
    plt.xticks(np.arange(40,101,20), fontsize = fs)
    plt.xlim(35,105)
plt.savefig(scra+variable+'_vert_seasons.jpg', bbox_inches='tight', dpi=300)
plt.savefig(med+variable+'_vert_seasons.jpg', bbox_inches='tight', dpi=300)
plt.savefig(med+variable+'_vert_seasons.pdf', bbox_inches='tight', dpi=300)
#plt.show()
plt.close()

