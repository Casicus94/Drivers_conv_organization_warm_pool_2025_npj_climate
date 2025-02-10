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
user = '/home/tompkins-archive2/acasallas/ERA5_real/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'

def hov_muller(dates,lat,lon,ntime):
    #Paths
    scr = '/home/netapp-clima/scratch/acasallas/'
    mdir = '/home/tompkins-archive/acasallas/ERA5/'
    user = '/home/tompkins-archive2/acasallas/ERA5_real/'
    hovchunk=81  # larger=smoother.
    #Constants
    Rd=287.05
    Cp=1005.0
    Lv=2.5e6
    Lf=334000

    #Start matrixes 
    sst_hov=np.zeros([1,ntime,nhov])
    lh_hov=np.zeros([1,ntime,nhov])
    sh_hov=np.zeros([1,ntime,nhov])
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
    tcwv_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/TCWV_area.nc')
    lh_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/LHF_area.nc')
    sh_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SHF_area.nc')
    lwtbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/LWTBNT_area.nc')
    lwctbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/LWCTBNT_area.nc')
    swtbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SWTBNT_area.nc')
    swctbn_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SWCTBNT_area.nc')
    ##LW
    lwt_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/TOA_net_LW_area.nc')
    lwct_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/TOA_net_LWC_area.nc')
    lwb_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SFC_net_LW_area.nc')
    lwcb_ds = xr.open_dataset(user+'Area_'+lat+'_'+lon+'/SFC_net_LWC_area.nc')

    #Starting calculation
    time1 = 0
    for itime,date in enumerate(dates):
        sst = np.array(sst_ds.sst.loc[date])
        tcwv = np.array(tcwv_ds.tcwv.loc[date])
        lh = np.array(lh_ds.slhf.loc[date])*(-1)/3600
        sh = np.array(sh_ds.sshf.loc[date])*(-1)/3600
        lwtbn = np.array(lwtbn_ds.LWTBNT.loc[date])
        lwctbn = np.array(lwctbn_ds.LWCTBNT.loc[date])
        swtbn = np.array(swtbn_ds.SWTBNT.loc[date])
        swctbn = np.array(swctbn_ds.SWCTBNT.loc[date])
        lwt = np.array(lwt_ds.ttr.loc[date])/3600
        lwct = np.array(lwct_ds.ttrc.loc[date])/3600
        lwb = np.array(lwb_ds.str.loc[date])/3600
        lwcb = np.array(lwcb_ds.strc.loc[date])/3600
        ###Sort the tcwv
        isort=np.argsort(tcwv.flatten()) # sort tcwv
        #Perturbations
        sstpert=sst.flatten()-np.mean(sst)
        lhpert=lh.flatten()-np.mean(lh)
        shpert=sh.flatten()-np.mean(sh)
        lwtbnpert = lwtbn.flatten()-np.mean(lwtbn)
        lwctbnpert = lwctbn.flatten()-np.mean(lwctbn)
        swtbnpert = swtbn.flatten()-np.mean(swtbn)
        swctbnpert = swctbn.flatten()-np.mean(swctbn)
        lwtpert = lwt.flatten()-np.mean(lwt)
        lwctpert = lwct.flatten()-np.mean(lwct)
        lwbpert = lwb.flatten()-np.mean(lwb)
        lwcbpert = lwcb.flatten()-np.mean(lwcb)
        #Fill the matrix
        sst_hov[0,itime,:]=np.mean(sstpert[isort].reshape(-1,hovchunk),axis=1)
        lh_hov[0,itime,:]=np.mean(lhpert[isort].reshape(-1,hovchunk),axis=1)
        sh_hov[0,itime,:]=np.mean(shpert[isort].reshape(-1,hovchunk),axis=1)
        lwtbn_hov[0,itime,:]=np.mean(lwtbnpert[isort].reshape(-1,hovchunk),axis=1)
        lwctbn_hov[0,itime,:]=np.mean(lwctbnpert[isort].reshape(-1,hovchunk),axis=1)
        swtbn_hov[0,itime,:]=np.mean(swtbnpert[isort].reshape(-1,hovchunk),axis=1)
        swctbn_hov[0,itime,:]=np.mean(swctbnpert[isort].reshape(-1,hovchunk),axis=1)
        tcwv_hov[0,itime,:]=np.mean(tcwv.flatten()[isort].reshape(-1,hovchunk),axis=1)
        lwt_hov[0,itime,:]=np.mean(lwtpert[isort].reshape(-1,hovchunk),axis=1)
        lwct_hov[0,itime,:]=np.mean(lwctpert[isort].reshape(-1,hovchunk),axis=1)
        lwb_hov[0,itime,:]=np.mean(lwbpert[isort].reshape(-1,hovchunk),axis=1)
        lwcb_hov[0,itime,:]=np.mean(lwcbpert[isort].reshape(-1,hovchunk),axis=1)
    return(sst_hov,tcwv_hov,lh_hov,sh_hov,lwtbn_hov,lwctbn_hov,swtbn_hov,swctbn_hov,lwt_hov,lwct_hov,lwb_hov,lwcb_hov)

####Dates
def sel_season(data,one,two,three,sta,med,end):
    tmp = data[data.index.month.isin([one,two,three,sta,med,end])]
    new_df = tmp.index.strftime('%Y-%m-%d %H:%M:%S').values
    return(new_df)

mdir1 = '/home/tompkins-archive/acasallas/'
latitude = '2-9'
longitude = '135-145'

dates_rev = pd.read_csv(mdir1+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/Dates_'+latitude+'_'+longitude+'_reversals.csv', parse_dates = True, index_col = 0)

dates_org = pd.read_csv(mdir1+'RRTMG_data/Area_'+latitude+'_'+longitude+'_ERA5/Dates_'+latitude+'_'+longitude+'_organised.csv', parse_dates = True, index_col = 0)

dat_rev = sel_season(dates_rev,12,1,2,3,4,5)
dat_org = sel_season(dates_org,12,1,2,3,4,5)

try:
    dates_tot = pd.read_csv(mdir1+'RRTMG_data/Dates.csv',parse_dates = True, index_col = 0)
except:
    dates_tot = [d.strftime('%Y-%m-%d %H:%M:%S') for d in pd.date_range('2016-01-01','2019-12-31', freq='1H')]
    dates_df = pd.DataFrame(dates_tot, columns=['Dates'])
    dates_df.to_csv(mdir1+'RRTMG_data/Dates.csv', index = False)
    dates_tot = pd.read_csv(mdir1+'RRTMG_data/Dates.csv',parse_dates = True, index_col = 0)

dat_mam = sel_season(dates_tot,12,1,2,3,4,5)
dat_jja = sel_season(dates_tot,6,7,8,9,10,11)

###We need to preprocess

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
ntime_rev = len(dat_rev)
ntime_org = len(dat_org)
ntime_mam = len(dat_mam)
ntime_jja = len(dat_jja)
#Location
lat,lon = '2-9','135-145'

try:
    var_rev = pd.read_csv(mdir1+'Data_to_plot/Hov_fluxes_rev_winter.csv')
    var_org = pd.read_csv(mdir1+'Data_to_plot/Hov_fluxes_org_winter.csv')
    var_jja = pd.read_csv(mdir1+'Data_to_plot/Hov_fluxes_summer.csv')
    var_mam = pd.read_csv(mdir1+'Data_to_plot/Hov_fluxes_winter.csv')
except:
    #Calculation
    print('Calculating Reversals')
    sst_rev,tcwv_rev,lh_rev,sh_rev,lwtbn_rev,lwctbn_rev,swtbn_rev,swctbn_rev,lwt_rev,lwct_rev,lwb_rev,lwcb_rev = hov_muller(dat_rev,lat,lon,ntime_rev)

    print('Calculating Organized')
    sst_org,tcwv_org,lh_org,sh_org,lwtbn_org,lwctbn_org,swtbn_org,swctbn_org,lwt_org,lwct_org,lwb_org,lwcb_org = hov_muller(dat_org,lat,lon,ntime_org)

    print('Calculating JJA')
    sst_jja,tcwv_jja,lh_jja,sh_jja,lwtbn_jja,lwctbn_jja,swtbn_jja,swctbn_jja,lwt_jja,lwct_jja,lwb_jja,lwcb_jja = hov_muller(dat_jja,lat,lon,ntime_jja)

    print('Calculating MAM')
    sst_mam,tcwv_mam,lh_mam,sh_mam,lwtbn_mam,lwctbn_mam,swtbn_mam,swctbn_mam,lwt_mam,lwct_mam,lwb_mam,lwcb_mam = hov_muller(dat_mam,lat,lon,ntime_mam)

    var_rev = [np.mean(sst_rev[0,:,:],axis=0),
               np.mean(lwctbn_rev[0,:,:],axis=0),np.mean(swctbn_rev[0,:,:],axis=0),
               np.mean(lh_rev[0,:,:],axis=0),np.mean(sh_rev[0,:,:],axis=0),
               np.mean(lwtbn_rev[0,:,:],axis=0)-np.mean(lwctbn_rev[0,:,:],axis=0),
               np.mean(swtbn_rev[0,:,:],axis=0)-np.mean(swctbn_rev[0,:,:],axis=0),
               np.mean(lh_rev[0,:,:],axis=0) + np.mean(sh_rev[0,:,:],axis=0) + np.mean(lwtbn_rev[0,:,:],axis=0) + np.mean(swtbn_rev[0,:,:],axis=0)]

    var_org = [np.mean(sst_org[0,:,:],axis=0),
               np.mean(lwctbn_org[0,:,:],axis=0),np.mean(swctbn_org[0,:,:],axis=0),
               np.mean(lh_org[0,:,:],axis=0),np.mean(sh_org[0,:,:],axis=0),
               np.mean(lwtbn_org[0,:,:],axis=0)-np.mean(lwctbn_org[0,:,:],axis=0),
               np.mean(swtbn_org[0,:,:],axis=0)-np.mean(swctbn_org[0,:,:],axis=0),
               np.mean(lh_org[0,:,:],axis=0) + np.mean(sh_org[0,:,:],axis=0) + np.mean(lwtbn_org[0,:,:],axis=0) + np.mean(swtbn_org[0,:,:],axis=0)]

    var_jja = [np.mean(sst_jja[0,:,:],axis=0), 
               np.mean(lwctbn_jja[0,:,:],axis=0),np.mean(swctbn_jja[0,:,:],axis=0),
               np.mean(lh_jja[0,:,:],axis=0),np.mean(sh_jja[0,:,:],axis=0),
               np.mean(lwtbn_jja[0,:,:],axis=0)-np.mean(lwctbn_jja[0,:,:],axis=0),
               np.mean(swtbn_jja[0,:,:],axis=0)-np.mean(swctbn_jja[0,:,:],axis=0),
               np.mean(lh_jja[0,:,:],axis=0) + np.mean(sh_jja[0,:,:],axis=0) + np.mean(lwtbn_jja[0,:,:],axis=0) + np.mean(swtbn_jja[0,:,:],axis=0)]

    var_mam = [np.mean(sst_mam[0,:,:],axis=0), 
               np.mean(lwctbn_mam[0,:,:],axis=0),np.mean(swctbn_mam[0,:,:],axis=0),
               np.mean(lh_mam[0,:,:],axis=0),np.mean(sh_mam[0,:,:],axis=0),
               np.mean(lwtbn_mam[0,:,:],axis=0)-np.mean(lwctbn_mam[0,:,:],axis=0),
               np.mean(swtbn_mam[0,:,:],axis=0)-np.mean(swctbn_mam[0,:,:],axis=0), 
               np.mean(lh_mam[0,:,:],axis=0) + np.mean(sh_mam[0,:,:],axis=0) + np.mean(lwtbn_mam[0,:,:],axis=0) + np.mean(swtbn_mam[0,:,:],axis=0)]
    ######## DataFrames
    test_rev = pd.DataFrame(np.column_stack((var_rev)), columns = ['SST','LW Clear','SW Clear','Latent Heat Flux','Sensible Heat Flux','LW Cloud','SW Cloud', 'Total'])
    test_rev.to_csv(mdir1+'Data_to_plot/Hov_fluxes_rev_winter.csv')
    ###
    test_org = pd.DataFrame(np.column_stack((var_org)), columns = ['SST','LW Clear','SW Clear','Latent Heat Flux','Sensible Heat Flux','LW Cloud','SW Cloud', 'Total'])
    test_org.to_csv(mdir1+'Data_to_plot/Hov_fluxes_org_winter.csv')
    ###
    test_mam = pd.DataFrame(np.column_stack((var_mam)), columns = ['SST','LW Clear','SW Clear','Latent Heat Flux','Sensible Heat Flux','LW Cloud','SW Cloud', 'Total'])
    test_mam.to_csv(mdir1+'Data_to_plot/Hov_fluxes_winter.csv')
    ###
    test_jja = pd.DataFrame(np.column_stack((var_jja)), columns = ['SST','LW Clear','SW Clear','Latent Heat Flux','Sensible Heat Flux','LW Cloud','SW Cloud', 'Total'])
    test_jja.to_csv(mdir1+'Data_to_plot/Hov_fluxes_summer.csv')
    print('Run me again, I created the files for plotting!')
    quit()

labels = ['Total','LW Clear','SW Clear','Latent Heat Flux','Sensible Heat Flux','LW Cloud','SW Cloud']
colors = ['k','steelblue','red','darkcyan','cyan','darkblue','darkred']

plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(2,2, left=0.05, right=0.95, hspace=0.2, wspace=0.15, top=0.9, bottom=0.1, width_ratios = [1,1])

ax = subplot(gs[0])
for i,var in enumerate(labels):
    plt.plot(hovx,var_jja[var], label = var, color = colors[i])
plt.xticks(np.arange(0,101,20))
plt.title('(a) Boreal Summer/Autumn',fontweight='bold', loc='left')
plt.yscale('symlog')
plt.ylabel('W m$^{-2}$')
plt.legend(frameon=False, ncol = 7, loc = 'upper center', bbox_to_anchor = (1.07,1.2))
plt.setp(ax.get_xticklabels(), visible=False)
plt.axhline(0,linestyle=':', color = 'grey')

ax = subplot(gs[1])
for i,var in enumerate(labels):
    plt.plot(hovx,var_mam[var], label = var, color = colors[i])
plt.xticks(np.arange(0,101,20))
plt.title('(b) Boreal Winter/Spring',fontweight='bold', loc='left')
plt.yscale('symlog')
plt.ylabel('W m$^{-2}$')
plt.setp(ax.get_xticklabels(), visible=False)
plt.axhline(0,linestyle=':', color = 'grey')

ax = subplot(gs[2])
for i,var in enumerate(labels):
    plt.plot(hovx,var_rev[var], label = var, color = colors[i])
plt.xticks(np.arange(0,101,20))
plt.title('(c) Boreal Winter/Spring Reversals',fontweight='bold', loc='left')
plt.yscale('symlog')
plt.ylabel('W m$^{-2}$')
plt.xlabel('TCWV %-tile')
plt.axhline(0,linestyle=':', color = 'grey')

ax = subplot(gs[3])
for i,var in enumerate(labels):
    plt.plot(hovx,var_org[var], label = var, color = colors[i])
plt.xticks(np.arange(0,101,20))
plt.title('(d) Boreal Winter/Spring Organized',fontweight='bold', loc='left')
plt.yscale('symlog')
plt.ylabel('W m$^{-2}$')
plt.xlabel('TCWV %-tile')
plt.axhline(0,linestyle=':', color = 'grey')

plt.savefig(scr+'hov_rev_org_fluxes_winter_summer.jpg', bbox_inches='tight', dpi=300)
plt.savefig(scr+'hov_rev_org_fluxes_winter_summer.pdf', bbox_inches='tight', dpi=300)
plt.show()
plt.close()

