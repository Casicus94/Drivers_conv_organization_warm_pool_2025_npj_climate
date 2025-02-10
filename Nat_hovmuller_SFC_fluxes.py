import xarray as xr
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.gridspec as gridspec
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from bisect import bisect_left
from scipy import stats
import Casicus as casi
import pandas as pd

scr = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/ERA5/'
user = '/home/netapp-clima/users/acasallas/ERA5_real/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'

def hov_muller(sta_date,end_date,lat,lon,ntime):
    tipo = 'Him'
    #Paths
    scr = '/home/netapp-clima/scratch/acasallas/'
    mdir = '/home/tompkins-archive/acasallas/ERA5/'
    user = '/home/netapp-clima/users/acasallas/ERA5_real/'
    obdr = '/home/tompkins-archive2/acasallas/Observations/' 
    hovchunk=81  #1024 # averaging length for chunking, power of 2, larger=smoother.
    #Constants
    Rd=287.05
    Cp=1005.0
    Lv=2.5e6
    Lf=334000

    #Start matrixes 
    times=np.array(range(ntime)) # start from day 1 for spin up
    sst_hov=np.zeros([1,ntime,nhov])
    lh_hov=np.zeros([1,ntime,nhov])
    tcwv_hov=np.zeros([1,ntime,nhov])
    lwtbn_hov = np.zeros([1,ntime,nhov])  
    lwctbn_hov = np.zeros([1,ntime,nhov])
    swtbn_hov = np.zeros([1,ntime,nhov])
    swctbn_hov = np.zeros([1,ntime,nhov])
    ##LW
    lwt_hov = np.zeros([1,ntime,nhov])
    lwct_hov = np.zeros([1,ntime,nhov])
    lwb_hov = np.zeros([1,ntime,nhov])
    lwcb_hov = np.zeros([1,ntime,nhov])
    #Read files
    sst_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SST_area.nc')
    sst_ds = sst_ds.sst.loc[sta_date:end_date]
    if tipo == 'ERA5':
        #sst_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SST_area.nc')
        #sst_ds = sst_ds.sst.loc[sta_date:end_date]
        tcwv_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/TCWV_area.nc')
        tcwv_ds = tcwv_ds.tcwv.loc[sta_date:end_date]
    else:
        #sst_ds = xr.open_dataset(obdr+'SST_Him_area_ERA5_grid.nc')
        #sst_ds = sst_ds.sea_surface_temperature.loc[sta_date:end_date]
        tcwv_ds = xr.open_dataset(obdr+'Mimic_area_ERA5_grid.nc')
        tcwv_ds = tcwv_ds.tpwGrid.loc[sta_date:end_date]
    lh_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/LHF_area.nc')
    lwtbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/LWTBNT_area.nc')
    lwctbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/LWCTBNT_area.nc')
    swtbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SWTBNT_area.nc')
    swctbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SWCTBNT_area.nc')
    ##LW
    lwt_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/TOA_net_LW_area.nc')
    lwct_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/TOA_net_LWC_area.nc')
    lwb_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SFC_net_LW_area.nc')
    lwcb_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SFC_net_LWC_area.nc')

    #Split dates
    lh_ds = lh_ds.slhf.loc[sta_date:end_date]
    lwtbn_ds = lwtbn_ds.LWTBNT.loc[sta_date:end_date] 
    lwctbn_ds = lwctbn_ds.LWCTBNT.loc[sta_date:end_date]
    swtbn_ds = swtbn_ds.SWTBNT.loc[sta_date:end_date]
    swctbn_ds =  swctbn_ds.SWCTBNT.loc[sta_date:end_date]
    lwt_ds = lwt_ds.ttr.loc[sta_date:end_date]
    lwct_ds = lwct_ds.ttrc.loc[sta_date:end_date]
    lwb_ds = lwb_ds.str.loc[sta_date:end_date]
    lwcb_ds = lwcb_ds.strc.loc[sta_date:end_date]
    #Starting calculation
    time1 = 0
    for itime in times:
        slice = itime + time1 
        tcwv = np.array(tcwv_ds[slice,:,:])
        sst = np.array(sst_ds[slice,:,:])
        lh = np.array(lh_ds[slice,:,:])*(-1)/3600
        lwtbn = np.array(lwtbn_ds[slice,:,:])
        lwctbn = np.array(lwctbn_ds[slice,:,:])
        swtbn = np.array(swtbn_ds[slice,:,:])
        swctbn = np.array(swctbn_ds[slice,:,:])
        lwt = np.array(lwt_ds[slice,:,:])/3600
        lwct = np.array(lwct_ds[slice,:,:])/3600
        lwb = np.array(lwb_ds[slice,:,:])/3600
        lwcb = np.array(lwcb_ds[slice,:,:])/3600
        ###Sort the tcwv
        isort=np.argsort(tcwv.flatten()) # sort tcwv
        #Perturbations
        sstpert=sst.flatten()-np.nanmean(sst)
        lhpert=lh.flatten()-np.mean(lh)
        lwtbnpert = lwtbn.flatten()-np.mean(lwtbn)
        lwctbnpert = lwctbn.flatten()-np.mean(lwctbn)
        swtbnpert = swtbn.flatten()-np.mean(swtbn)
        swctbnpert = swctbn.flatten()-np.mean(swctbn)
        lwtpert = lwt.flatten()-np.mean(lwt)
        lwctpert = lwct.flatten()-np.mean(lwct)
        lwbpert = lwb.flatten()-np.mean(lwb)
        lwcbpert = lwcb.flatten()-np.mean(lwcb)
        #Fill the matrix
        sst_hov[0,itime,:]=np.nanmean(sstpert[isort].reshape(-1,hovchunk),axis=1)
        lh_hov[0,itime,:]=np.mean(lhpert[isort].reshape(-1,hovchunk),axis=1)
        lwtbn_hov[0,itime,:]=np.mean(lwtbnpert[isort].reshape(-1,hovchunk),axis=1)
        lwctbn_hov[0,itime,:]=np.mean(lwctbnpert[isort].reshape(-1,hovchunk),axis=1)
        swtbn_hov[0,itime,:]=np.mean(swtbnpert[isort].reshape(-1,hovchunk),axis=1)
        swctbn_hov[0,itime,:]=np.mean(swctbnpert[isort].reshape(-1,hovchunk),axis=1)
        tcwv_hov[0,itime,:]=np.mean(tcwv.flatten()[isort].reshape(-1,hovchunk),axis=1)
        lwt_hov[0,itime,:]=np.mean(lwtpert[isort].reshape(-1,hovchunk),axis=1)
        lwct_hov[0,itime,:]=np.mean(lwctpert[isort].reshape(-1,hovchunk),axis=1)
        lwb_hov[0,itime,:]=np.mean(lwbpert[isort].reshape(-1,hovchunk),axis=1)
        lwcb_hov[0,itime,:]=np.mean(lwcbpert[isort].reshape(-1,hovchunk),axis=1)
    return(sst_hov,tcwv_hov,lh_hov,lwtbn_hov,lwctbn_hov,swtbn_hov,swctbn_hov,lwt_hov,lwct_hov,lwb_hov,lwcb_hov)

###First we need to preprocess

print('Starting Hovmuller script!')
qv_press=1100. # low press qv integration hPA

cmap="bwr"
cmap_r="bwr_r"

ylab='Time (days)'

xlab= "TCWV %-tile"
funits="(W m$^{-2}$)"

sstconts=[-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0]
sstconts=[-0.2,-0.1,-0.05,0.05,0.1,0.2]

#Domain and splits
nx=57
ny=81

npts=nx*ny
hovchunk=81  #1024 # averaging length for chunking, power of 2, larger=smoother.
nhov=int(npts/hovchunk)
hovx=100*(np.arange(nhov)+0.5)/nhov
cl=[0.0]
fs=9
vmax = 0.5

#Dates 
sta_date,end_date = '2017-03-15','2017-05-15'
#sta_date,end_date = '2017-06-15','2017-08-15'
tipo = 'Him'
time_tot = [d.strftime('%Y-%m-%d') for d in pd.date_range(sta_date,end_date, freq='1d')]
dates_plot = [d.strftime('%Y-%m-%d') for d in pd.date_range(sta_date,end_date, freq='3d')]
tval_plot = pd.date_range(sta_date,end_date+' 23:00:00', freq='1H')
ntime = len(time_tot)*24
#Location
lat,lon = '2-9','135-145'

#Calculation
sst_hov,tcwv_hov,lh_hov,lwtbn_hov,lwctbn_hov,swtbn_hov,swctbn_hov,lwt_hov,lwct_hov,lwb_hov,lwcb_hov = hov_muller(sta_date,end_date,lat,lon,ntime)

if tipo == 'Him':
    sst_hov = sst_hov*3
    levels = [np.arange(-1,1.01,0.25),np.arange(-30,31,5),np.arange(-30,31,5),
              np.arange(-30,31,5),np.arange(-30,31,5),np.arange(-30,31,5)]
else:
    levels = [np.arange(-0.4,0.41,0.1),np.arange(-30,31,5),np.arange(-30,31,5),
              np.arange(-30,31,5),np.arange(-30,31,5),np.arange(-30,31,5)]

varis = [sst_hov,lwctbn_hov,swctbn_hov,lh_hov,lwtbn_hov-lwctbn_hov,swtbn_hov-swctbn_hov]
#titles = ['(g) SST','(h) LW Clear','(i) SW Clear','(j) Latent Heat Flux','(k) LW Cloud','(l) SW Cloud']
titles = ['(a) SST','(b) LW Clear','(c) SW Clear','(d) Latent Heat Flux','(e) LW Cloud','(f) SW Cloud']

plt.figure(figsize=(14,16))
gs = gridspec.GridSpec(2,3, left=0.13, right=0.95, hspace=0.15, wspace=0.1, top=0.95, bottom=0.1)

pos = [0,1,2,3,4,5]
for i,var in enumerate(varis):
    ax = subplot(gs[pos[i]])
    im = plt.contourf(hovx,tval_plot, var[0], cmap=cmap, levels=levels[i], extend='both')
    plt.colorbar(im)
    ax.set_yticks(dates_plot)
    ax.set_yticklabels(dates_plot)
    plt.xlabel('TCWV %-tile')
    plt.title(titles[i],fontweight='bold', loc='left')
    if i == 1 or i == 2 or i == 4 or i == 5:
        plt.setp(ax.get_yticklabels(), visible=False)

#ax2 = subplot(gs[11])
#plt.colorbar(im, cax=ax2)

#ax1 = subplot(gs[0,1:3])
#cbar1=plt.colorbar(im, cax=ax1, orientation='horizontal')
#cbar1.ax.xaxis.set_ticks_position("top")

plt.savefig(med+'hov_tot_fluxes_'+sta_date+'_'+end_date+'_'+tipo+'.jpg', bbox_inches='tight', dpi=300)
plt.savefig(med+'hov_tot_fluxes_'+sta_date+'_'+end_date+'_'+tipo+'.pdf', bbox_inches='tight', dpi=300)
plt.show()
plt.close()

#######################################################################
#                             LW TOP and SFC                          #
#######################################################################
varis = [lwct_hov,lwcb_hov,lwt_hov-lwct_hov,lwb_hov-lwcb_hov]
levels = [np.arange(-30,31,5),np.arange(-30,31,5),
          np.arange(-30,31,5),np.arange(-30,31,5)]
titles = ['(a) LW TOP Clear','(b) LW SFC Clear',
          '(c) LW TOP Cloud','(d) LW SFC Cloud',]

plt.figure(figsize=(9,18))
gs = gridspec.GridSpec(2,2, left=0.13, right=0.95, hspace=0.25, wspace=0.1, top=0.9, bottom=0.1)

for i,var in enumerate(varis):
    ax = subplot(gs[i])
    im = plt.contourf(hovx,tval_plot, var[0], cmap=cmap, levels=levels[i], extend='both')
    plt.colorbar(im)
    ax.set_yticks(dates_plot)
    ax.set_yticklabels(dates_plot)
    plt.xlabel('TCWV %-tile')
    plt.title(titles[i],fontweight='bold', loc='left')
    if i == 1 or i == 3 or i == 4:
        plt.setp(ax.get_yticklabels(), visible=False)

plt.savefig(med+'hov_LW_tot_'+sta_date+'_'+end_date+'_'+tipo+'.jpg', bbox_inches='tight', dpi=300)
plt.savefig(med+'hov_LW_tot_'+sta_date+'_'+end_date+'_'+tipo+'.pdf', bbox_inches='tight', dpi=300)
#plt.show()


