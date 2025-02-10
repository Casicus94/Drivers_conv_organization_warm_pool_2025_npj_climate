from netCDF4 import Dataset
import pdb
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *

def hov_fluxes(name):
    flxlist=[
    # Fluxes
    {"nm":"ACHFX","title":"sensible heat flux"},
    {"nm":"ACLHF","title":"latent heat flux"},

    # Clear sky
    {"nm":"ACLWBNETC","title":"LW Clear net surface"},
    {"nm":"ACLWTBNETC","title":"LW Clear net atmos"},
    {"nm":"ACLWTNETC","title":"LW Clear net TOA"},
    {"nm":"ACSWBNETC","title":"SW Clear net surface"},
    {"nm":"ACSWTBNETC","title":"SW Clear net atmosphere"},
    {"nm":"ACSWTNETC","title":"SW Clear net TOA"},

    # Net full 
    {"nm":"ACLWBNET","title":"LW net surface"},
    {"nm":"ACLWTBNET","title":"LW net atmos"},
    {"nm":"ACLWTNET","title":"LW net TOA"},
    {"nm":"ACSWBNET","title":"SW net surface"},
    {"nm":"ACSWTBNET","title":"SW net atmosphere"},
    {"nm":"ACSWTNET","title":"SW net TOA"},
    ]

    if name == 'JJA' or name == 'u_20mar' or name == 'v_20mar':
        tcwv_ds = xr.open_dataset(obdr+name+'/qv_vert_'+name+'.nc')
        fluxfile = obdr+name+'/allfld_test_'+name+'_flux_d3600.nc'
        sst_ds = xr.open_dataset(obdr+name+'/TSK_'+name+'.nc')
    else:
        tcwv_ds = xr.open_dataset(reals+name+'/qv_vert_'+name+'.nc')
        fluxfile = reals+name+'/allfld_test_'+name+'_flux_d3600.nc'
        sst_ds = xr.open_dataset(reals+name+'/TSK_'+name+'.nc')
    time1 = 0*24
    time2 = len(tcwv_ds.QVAPOR) - 1
    ntime = time2-time1+1
    dtime = 3600 #This depends on the ouput frequency
    difftime = str(dtime)
    fstr = difftime+'-'+str(time1)+'-'+str(time2)
    nx = 552; ny = 399; nz = 32; npts = nx*ny # hardwired for moment
    hovchunk=1656  #1024 # averaging length for chunking, power of 2, larger=smoother.
    nhov=int(npts/hovchunk)
    hovx=100*(np.arange(nhov)+0.5)/nhov
    ilev1=0
    ilev2=8
    #constants
    Rd = 287.05
    Cp = 1005.0
    Lv = 2.5e6
    Cpdsinday = Cp/86400.
    #Start matrixes 
    times = np.array(range(ntime)) # start from day 1 for spin up
    hovtimes = (times+time1)*dtime/86400.
    tcwv_hov = np.zeros([1,ntime,nhov])
    sst_hov = np.zeros([1,ntime,nhov])
    lh_hov = np.zeros([1,ntime,nhov])
    hov={}
    for iflx in flxlist:
        hov[iflx["nm"]] = np.zeros([1,ntime,nhov])

    ds0=Dataset(fluxfile)

    for itime in times: 
        print('Time: '+str(itime))
        lh = np.array(ds0.variables["ACLHF"][itime,:,:])
        tcwv = np.array(tcwv_ds.variables['QVAPOR'][itime,:,:])
        sst =  np.array(sst_ds.variables['TSK'][itime,:,:])
        isort = np.argsort(tcwv.flatten())
        sstpert = sst.flatten()-np.mean(sst)
        sst_hov[0,itime,:] = np.mean(sstpert[isort].reshape(-1,hovchunk),axis=1)
        tcwv_hov[0,itime,:] = np.mean(tcwv.flatten()[isort].reshape(-1,hovchunk),axis=1)
    
        for iflx in flxlist:
            flxname=iflx["nm"]
            flx=np.array(ds0.variables[flxname][itime,:,:]).flatten()
            flx-=flx.mean() # anomaly
            # sort according to TCWV and then chunk-mean
            hov1=np.mean(flx[isort].reshape(-1,hovchunk),axis=1)
            # store to dictory slice
            hov[flxname][0,itime,:]=hov1
    return(hov,hovx,ntime)

### Paths!
scra = '/home/netapp-clima/scratch/acasallas/'
reals = '/home/tompkins-archive/acasallas/Real1_run/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'
mdir = '/home/tompkins-archive/acasallas/'
obdr = '/home/tompkins-archive2/acasallas/Real1_run/'

names = ['control','JJA','20_mar','qv_20mar','uv_20mar','u_20mar','v_20mar','12_apr','qv_12apr','uv_12apr','u_12apr','v_12apr']
data = {}

for name in names:
    print('---------> Running: '+name)
    hov, hovx, ntime = hov_fluxes(name) 
    data[name] = hov
    if name == 'JJA':
        dates = [d.strftime('%Y-%m-%d') for d in pd.date_range('2017-06-15 00:00:00','2017-08-15 00:00:00', freq='3d')]
        tval = pd.date_range('2017-06-15 00:00:00','2017-08-15 00:00:00', freq='1H')
    else:
        dates = [d.strftime('%Y-%m-%d') for d in pd.date_range('2017-03-15 00:00:00','2017-05-13 23:00:00', freq='3d')]
        tval = pd.date_range('2017-03-15 00:00:00','2017-05-13 23:00:00', freq='1H')
    varis = ['ACLHF','ACLWTBNETC','ACSWTBNETC','ACHFX','ACLWTBNET','ACSWTBNET']
    title = ['(a) Latent Heat Flux','(b) LW Clear','(c) SW Clear','(d) Sensible Heat Flux','(e) LW Cloud','(f) SW Cloud']
    levs = np.arange(-30,30.1,5)
    fig = plt.figure(figsize=(10,10))
    gs = GridSpec(2,4,left = 0.1, right = 0.92, hspace=0.15, wspace=0.15, top = 0.9, bottom = 0.08, height_ratios = [1,1], width_ratios = [1,1,1,0.1])
    for i,var in enumerate(varis):
        if i < 3:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i+1])
        if var == 'ACLWTBNET' or var == 'ACSWTBNET':
            im = plt.contourf(hovx,tval, hov['AC'+var[2:4]+'TBNET'][0,:len(tval),:]-hov['AC'+var[2:4]+'TBNETC'][0,:len(tval),:], cmap = 'bwr', levels = levs, extend = 'both')
        else:
            im = plt.contourf(hovx,tval, hov[var][0,:len(tval),:], cmap = 'bwr', levels = levs, extend = 'both')
        ax.set_yticks(dates)
        ax.set_yticklabels(dates)
        plt.xticks(np.arange(0,100.1,20))
        plt.title(title[i], loc = 'left')
        plt.ylabel('W m$^{-2}$')
        if i > 2:
            plt.xlabel('TCWV %-tile')
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        if i == 1 or i == 2 or i == 4 or i == 5:
            plt.setp(ax.get_yticklabels(), visible=False)
    ax = plt.subplot(gs[:,3])
    cbar = plt.colorbar(im, cax = ax, orientation = 'vertical', shrink = 0.4, ticks = np.arange(-30,30.1,10))
    #plt.savefig(med+name+'/Rad_SF_Fluxes_'+name+'.jpg', bbox_inches = 'tight', dpi = 300)
    #plt.savefig(med+name+'/Rad_SF_Fluxes_'+name+'.pdf', bbox_inches = 'tight', dpi = 300)
    #plt.show()
    plt.close('all')

### Means
means = {}
varis = ['ACLHF','ACLWTBNETC','ACSWTBNETC','ACHFX','ACLWTBNET','ACSWTBNET']

for name in names:
    for i,var in enumerate(varis):
        if i == 0:
            means[name] = data[name][var]
        else:
            means[name] = means[name] + data[name][var] 

col1 = ['darkcyan','red','royalblue','plum','mediumorchid','blueviolet','purple']
col2 = ['darkcyan','red','royalblue','plum','mediumorchid','blueviolet','purple']
lab1 = ['MAM','JJA','Adv. Moist Org.','Moist Org.','Wind Org.','Zonal Org.','Meridional Org.']
lab2 = ['MAM','JJA','Adv. Moist Rev.','Moist Rev.','Wind Rev.','Zonal Rev.','Meridional Rev.']

name1 = ['control','JJA','20_mar','qv_20mar','uv_20mar','u_20mar','v_20mar']
name2 = ['control','JJA','12_apr','qv_12apr','uv_12apr','u_12apr','v_12apr']

fig = plt.figure(figsize=(8,6))
gs = GridSpec(2,1,left = 0.1, right = 0.92, hspace=0.15, wspace=0.15, top = 0.9, bottom = 0.08)
ax = plt.subplot(gs[0])
for i,name in enumerate(name1):
    plt.plot(np.linspace(0,100,len(np.mean(means[name][0],axis=0))),np.mean(means[name][0],axis=0), color = col1[i], label = lab1[i])
    plt.title('(a)', loc = 'left')
    plt.axhline(0, linestyle = ':', color = 'k')
    plt.legend(frameon = False, ncol = 3)
    ax.xaxis.set_major_formatter(NullFormatter())    
    plt.ylabel('W m$^{-2}$')

ax = plt.subplot(gs[1])
for i,name in enumerate(name2):
    plt.plot(np.linspace(0,100,len(np.mean(means[name][0],axis=0))),np.mean(means[name][0],axis=0), color = col2[i], label = lab2[i])
    plt.title('(b)', loc = 'left')
    plt.axhline(0, linestyle = ':', color = 'k')
    plt.legend(frameon = False, ncol = 3)
plt.xlabel('TCWV %-tile')
plt.ylabel('W m$^{-2}$')
plt.savefig(scra+'/Total_Fluxes.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(scra+'/Total_Fluxes.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()





