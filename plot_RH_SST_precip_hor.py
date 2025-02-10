import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *
import pdb

### Directories
mdir = '/home/tompkins-archive/acasallas/'
user = '/home/tompkins-archive2/acasallas/ERA5_real/'
scra =  '/home/netapp-clima/scratch/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'
obdr = '/home/tompkins-archive2/acasallas/Observations/'

###Data
tipo = 'ERA5' # ERA5 or Himawari

if tipo == 'ERA5':
    ds_sst = xr.open_dataset(user+'SST_area_zonmean.nc')
    wei_sst = np.cos(np.deg2rad(ds_sst.sst.lat)) / np.cos(np.deg2rad(ds_sst.sst.lat)).mean(dim='lat')
    sst_sea = (ds_sst.sst).groupby('time.season').mean()
elif tipo == 'Himawari':
    ds_sst = xr.open_dataset(obdr+'SST_Him_area_ERA5_grid.nc')
    ds_sst = ds_sst.mean(dim=['longitude'])
    wei_sst = np.cos(np.deg2rad(ds_sst.sea_surface_temperature.latitude)) / np.cos(np.deg2rad(ds_sst.sea_surface_temperature.latitude)).mean(dim='latitude')
    sst_sea = (ds_sst.sea_surface_temperature).groupby('time.season').mean()

ds_pr = xr.open_dataset(user+'Precip_area_zonmean.nc')
ds_rh = xr.open_dataset(user+'rela_humidity_area_zonmean.nc')

##Seasons
ssts=[sst_sea[0],sst_sea[2],sst_sea[1],sst_sea[3]]

wei_pr = np.cos(np.deg2rad(ds_pr.Precip.lat)) / np.cos(np.deg2rad(ds_pr.Precip.lat)).mean(dim='lat')
pr_sea = (ds_pr.Precip).groupby('time.season').sum()/6
prs=[pr_sea[0],pr_sea[2],pr_sea[1],pr_sea[3]]

wei_rh = np.cos(np.deg2rad(ds_rh.r.lat)) / np.cos(np.deg2rad(ds_rh.r.lat)).mean(dim='lat')
rh_sea = (ds_rh.r).groupby('time.season').mean()
rhs=[rh_sea[0],rh_sea[2],rh_sea[1],rh_sea[3]]

##### Starting plot
seasons = ['(e) DJF','(f) MAM','(g) JJA','(h) SON']
x = np.arange(2,9.01,0.1)
xs = [0.15,0.17,0.13,0.16]
titulos = ['(a) DJF','(b) MAM','(c) JJA','(d) SON']

fig = plt.figure(figsize=(16,6)) 
gs = GridSpec(2,5,left = 0.09, right = 0.95, hspace=0.1, wspace=0.85, top = 0.95, bottom = 0.08, height_ratios = [1,0.75], width_ratios = [1.2,1.2,1.2,1.2,0.1])

for i,rh in enumerate(rhs):
    ax = plt.subplot(gs[0,i])
    im = plt.contourf(rh.lat,rh.level,rh.mean(dim='lon'), cmap = 'YlGnBu', levels=np.arange(40,100.01,2.5), extend = 'both')
    plt.ylabel('Pressure (hPa)')
    #plt.xlabel('Latitude')
    plt.ylim(1000,100)
    plt.xticks(np.arange(2,9.1,1))
    ax.xaxis.set_major_formatter(NullFormatter())
    plt.title(titulos[i], fontsize =10, loc = 'left')

gs1 = GridSpec(1,5)
gs1.update(left=0.8, right = 0.9, top= 0.95, bottom=0.475)
ax = plt.subplot(gs1[0,4])
cbar = plt.colorbar(im, cax = ax, ticks = np.arange(40,100.01,5))
cbar.set_label('Relative Humidity (%)') 


for i,sst in enumerate(ssts):
    ax = plt.subplot(gs[1,i])
    if tipo == 'ERA5':
        ax.plot(sst.lat.values,sst, color = 'darkred')
        plt.ylim(301.75,303.5)
    else:
        ax.plot(sst.latitude.values,sst, color = 'darkred')
        plt.ylim(301.75,303.85)
    if i >= 0:
        ax.set_ylabel('SST (K)', color = 'darkred')
    else:
        ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.label.set_color('darkred')
    plt.xticks(np.arange(2,9.1,1))
    plt.xlabel('Latitude')
    plt.xlim(2,9)
    ax.tick_params(axis='y', colors = 'darkred')
    ax2 = ax.twinx()
    ax2.plot(pr_sea.lat.values,prs[i], color = 'darkblue')
    plt.yticks(np.arange(400,1200.01,200))
    plt.ylim(300,1320)
    #plt.ylim(20,60)
    #plt.yticks(np.arange(25,55.01,5))
    plt.xticks(np.arange(2,9.1,1))
    plt.xlabel('Latitude')
    ax2.set_ylabel('Precip. (mm)', color = 'darkblue')
    ax2.yaxis.label.set_color('darkblue')
    ax2.tick_params(axis='y', colors = 'darkblue')
    plt.title(seasons[i], loc= 'left', y=0.88, fontsize = 10)
#plt.savefig(med+'SST_hum_rain_sea_hor_'+tipo+'.jpg', bbox_inches='tight', dpi = 300)
#plt.savefig(med+'SST_hum_rain_sea_hor_'+tipo+'.pdf', bbox_inches='tight', dpi = 300)
#plt.savefig(scra+'SST_hum_rain_sea_hor_'+tipo+'.jpg', bbox_inches='tight', dpi = 300)
#plt.savefig(scra+'SST_hum_rain_sea_hor_'+tipo+'.pdf', bbox_inches='tight', dpi = 300)
#plt.show()
plt.close()

###################################
##### Organized and reversals #####
###################################
print('Organized and Reversal cases!')
def sum_them(ds_tcwv,ds_sst,ds_rh,dates):
    datasets = {'Precip':[], 'SST':[], 'RH':[]}
    for i,date in enumerate(dates):
        datasets['Precip'].append(ds_tcwv.Precip.loc[date[0:10]])
        if tipo == 'ERA5': 
            datasets['SST'].append(ds_sst.sst.loc[date])
        else:
            datasets['SST'].append(ds_sst.sea_surface_temperature.loc[date])
        datasets['RH'].append(ds_rh.r.loc[date])
    ###Means
    #Precip
    pr_rev = xr.concat(datasets['Precip'], dim='time')
    #SST
    sst_rev = xr.concat(datasets['SST'], dim='time')
    #RH
    rh_rev = xr.concat(datasets['RH'], dim='time')
    return(pr_rev,sst_rev,rh_rev)

def sel_season(data,sta,med,end):
    tmp = data[data.index.month.isin([sta,med,end])]
    new_df = tmp.index.strftime('%Y-%m-%d %H:%M:%S').values
    return(new_df)

def sel_boreal(data,dic,ene,feb,ini,med,end):
    tmp = data[data.index.month.isin([dic,ene,feb,ini,med,end])]
    new_df = tmp.index.strftime('%Y-%m-%d %H:%M:%S').values
    return(new_df)

### Dates
import pandas as pd
import sys
latitude = '2-9'
longitude = '135-145'

dates_rev = pd.read_csv(mdir+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/Dates_'+latitude+'_'+longitude+'_reversals.csv', parse_dates = True, index_col = 0)
dates_org = pd.read_csv(mdir+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/Dates_'+latitude+'_'+longitude+'_organised.csv', parse_dates = True, index_col = 0)

dat_rev = dates_rev.index.strftime('%Y-%m-%d %H:%M:%S').values
dat_org = dates_org.index.strftime('%Y-%m-%d %H:%M:%S').values

ini = sys.argv[1]
med = sys.argv[2]
fin = sys.argv[3]
dat_org = sel_boreal(dates_org,12,1,2,3,4,5)
dat_rev = sel_boreal(dates_rev,12,1,2,3,4,5)
#dat_org = sel_season(dates_org,1,2,3)
#dat_rev = sel_season(dates_rev,1,2,3)

print('##### Selecting reversal days #####')
pr_rev,sst_rev,rh_rev = sum_them(ds_pr,ds_sst,ds_rh,dat_rev)
print('##### Selecting organized days #####')
pr_org,sst_org,rh_org = sum_them(ds_pr,ds_sst,ds_rh,dat_org)

prs = [pr_org.mean(dim=['time','lon']),pr_rev.mean(dim=['time','lon'])]
if tipo == 'ERA5':
    ssts = [sst_org.mean(dim = ['time','lon']),sst_rev.mean(dim = ['time','lon'])] 
else:
    ssts = [sst_org.mean(dim = ['time']),sst_rev.mean(dim = ['time'])]
rhs = [rh_org.mean(dim=['time','lon']),rh_rev.mean(dim=['time','lon'])]
title = ['(c) Org Boreal Winter/Spring','(d) Rev Boreal Winter/Spring']
tituls = ['(a) Org Boreal Winter/Spring','(b) Rev Boreal Winter/Spring']
#title = ['(c) Org DJF','(d) Rev DJF']
#tituls = ['(a) Org DJF','(b) Rev DJF']

fig = plt.figure(figsize=(8,6)) #8,16
gs = GridSpec(2,3,left = 0.1, right = 0.95, hspace=0.1, wspace=0.875, top = 0.95, bottom = 0.08, height_ratios = [1,0.75], width_ratios = [1,1,0.2])

for i,rh in enumerate(rhs):
    ax = plt.subplot(gs[0,i])
    im = plt.contourf(rh.lat,rh.level,rh, cmap = 'YlGnBu', levels=np.arange(40,100.01,2.5), extend = 'both')
    plt.ylabel('Pressure (hPa)')
    plt.xticks(np.arange(2,9.01,1))
    plt.ylim(1000,100)
    plt.title(tituls[i], fontsize =9, loc = 'left')
    ax.xaxis.set_major_formatter(NullFormatter())

gs1 = GridSpec(1,3)
gs1.update(left=0.7, right = 0.8, top= 0.95, bottom=0.475)
ax = plt.subplot(gs1[0,2])
cbar = plt.colorbar(im, cax=ax, ticks = np.arange(40,100.01,5))
cbar.set_label('Relative Humidity (%)')

for i,sst in enumerate(ssts):
    ax = plt.subplot(gs[1,i])
    if tipo == 'ERA5':
        ax.plot(sst.lat.values,sst, color = 'darkred')
        plt.ylim(301.8,303.15)
        plt.yticks(np.arange(302.0,303.01,0.25))
    else:
        ax.plot(sst.latitude.values,sst, color = 'darkred')
        plt.ylim(301.75,303.5)
        plt.yticks(np.arange(301.75,303.51,0.25))
    ax.set_ylabel('SST (K)', color = 'darkred')
    ax.yaxis.label.set_color('darkred')
    plt.xticks(np.arange(2,9.1,1))
    plt.xlabel('Latitude')
    plt.xlim(2,9)
    ax.tick_params(axis='y', colors = 'darkred')
    ax2 = ax.twinx()
    ax2.plot(pr_sea.lat.values,prs[i], color = 'darkblue')
    #plt.yticks(np.arange(400,1200.01,200))
    #plt.ylim(2.5,67.5)
    plt.ylim(10,60)
    #plt.yticks(np.arange(10,60.01,10))
    plt.yticks(np.arange(10,60.01,10))
    plt.xticks(np.arange(2,9.1,1))
    plt.xlabel('Latitude')
    ax2.set_ylabel('Precip. (mm)', color = 'darkblue')
    ax2.yaxis.label.set_color('darkblue')
    ax2.tick_params(axis='y', colors = 'darkblue')
    plt.title(title[i], loc='left', y=0.88, fontsize = 9)
#plt.savefig(med+'SST_rain_RH_rev_org_MAM_hor_'+tipo+'.jpg', bbox_inches='tight', dpi = 300)
#plt.savefig(med+'SST_rain_RH_rev_org_MAM_hor_'+tipo+'.pdf', bbox_inches='tight', dpi = 300)
plt.savefig(scra+'SST_rain_RH_rev_org_MAM_hor_'+tipo+'.jpg', bbox_inches='tight', dpi = 300)
plt.savefig(scra+'SST_rain_RH_rev_org_MAM_hor_'+tipo+'.pdf', bbox_inches='tight', dpi = 300)

plt.show()
plt.close()

