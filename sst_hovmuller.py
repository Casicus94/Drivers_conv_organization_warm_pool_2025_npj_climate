from netCDF4 import Dataset
import pdb
import xarray as xr
import glob, os
import matplotlib as mpl
import numpy as np
import pandas as pd
from bisect import bisect_left
from joblib import dump,load
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import Casicus as casi

#Paths!
scra = '/home/netapp-clima/scratch/acasallas/wrf/Real_run/WRF/'
mdir = '/home/tompkins-archive/acasallas/'
reals = '/home/tompkins-archive/acasallas/Real1_run/'
user = '/home/netapp-clima/users/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'

print('Starting Hovmuller script!')
time1 = 0*24
time2 = 55*24
ntime=time2-time1+1 # 240 quick hack should be read from the netcdf file really...
yticki= 5.0
qv_press=1100. # low press qv integration hPA

cvar = 'v_12apr'

runlist=['test_'+cvar+'.nc']

nrun=len(runlist)

cmap="bwr"
cmap_r="bwr_r"

dtime=3600 #This depends on the ouput frequency
difftime=str(dtime)

ylab='Time (days)'

xlab= "TCWV %-tile"
funits="(W m$^{-2}$)"
fstr=difftime+'-'+str(time1)+'-'+str(time2)

sstconts=[-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0]
sstconts=[-0.2,-0.1,-0.05,0.05,0.1,0.2]

# fluxlist is a list of dictionaries
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

#Domain and splits
nx=552
ny=399
nz=32
npts=nx*ny # hardwired for moment
hovchunk=1656  #1024 # averaging length for chunking, power of 2, larger=smoother.
nhov=int(npts/hovchunk)
hovx=100*(np.arange(nhov)+0.5)/nhov
ilev1=0
ilev2=8

#constants
Rd=287.05
Cp=1005.0
Lv=2.5e6
Cpdsinday=Cp/86400.

#Start matrixes 
times=np.array(range(ntime)) # start from day 1 for spin up
#slices=times+time1
hovtimes=(times+time1)*dtime/86400.

heights=np.array(range(nz))
temp_hov=np.zeros([nrun,nz,ntime,nhov])
qv_hov=np.zeros([nrun,nz,ntime,nhov])
qi_hov=np.zeros([nrun,nz,ntime,nhov])
ql_hov=np.zeros([nrun,nz,ntime,nhov])
wp_hov=np.zeros([nrun,nz,ntime,nhov])
sf_hov=np.zeros([nrun,nz,ntime,nhov])
lw_hovf=np.zeros([nrun,nz,ntime,nhov])
sw_hovf=np.zeros([nrun,nz,ntime,nhov])
tcwv_hov=np.zeros([nrun,ntime,nhov])
sst_hov=np.zeros([nrun,ntime,nhov])
lh_hov=np.zeros([nrun,ntime,nhov])
lts_hov=np.zeros([nrun,ntime,nhov])

print("Read ",ntime,"timeslices")
hov={}
for iflx in flxlist:
    hov[iflx["nm"]]=np.zeros([nrun,ntime,nhov])

for irun,run in enumerate(runlist):
    difftimefstr=fstr+"_tcwv"+str(qv_press)
    #reset ds?
    ds0=ds1=None
    dfile=reals+cvar+'/'+run
    #dfile=scra+run
    print("opening ",dfile)
    ds1=Dataset(dfile)
    #fluxfile=mdir+run+'_run/allfld_'+run+'_flux_d'+difftime+'_verified.nc'
    fluxfile=reals+cvar+'/allfld_test_'+cvar+'_flux_d3600.nc'
    #fluxfile=scra+'allfld_test_'+cvar+'_flux_d3600.nc'
    ds0=Dataset(fluxfile)
    print ("starting step loop")
    for itime in times:
        itime1 = itime+time1
        slice=itime1+time1
        print("slice ",itime1,slice)
        print("reading 3D vars")
        # 3D
        qv = np.array(ds1.variables["QVAPOR"][slice,:,:,:])
        ql = np.array(ds1.variables["QCLOUD"][slice,:,:,:])
        qi = np.array(ds1.variables["QICE"][slice,:,:,:])
        wvel = np.array(ds1.variables["W"][slice,:,:,:])
        theta = np.array(ds1.variables["T"][slice,:,:,:])
        pr = np.array(ds1.variables["P"][slice,:,:,:])
        pb = np.array(ds1.variables["PB"][slice,:,:,:])
        ph = np.array(ds1.variables["PH"][slice,:,:,:])
        phb = np.array(ds1.variables["PHB"][slice,:,:,:])

        print("reading 2D vars")

        # 2D
        ps=   np.array(ds1.variables["PSFC"][slice,:,:])
        sst=  np.array(ds1.variables["TSK"][slice,:,:])
        lh=   np.array(ds0.variables["ACLHF"][slice,:,:])
        #pdb.set_trace()
        znu=  np.array(ds1.variables["ZNU"][0,:])*ps[0,0]/100.
        lw_hgt=ds0.variables["ACLWTNET"][slice,:,:]
        sw_hgt=ds0.variables["ACSWTNET"][slice,:,:]
        press=pr[:,:,:]+pb
        hgt=(ph[:,:,:]+phb[:,:,:])/9.81
        theta = theta + 300
        temp = theta*np.power(press/1.e5,Rd/Cp)
        tempv=(1.0+0.61*qv)*temp # tempv is vtemp

        #Averages
        avpress=np.mean(press,axis=(1,2))
        avtemp=np.mean(temp,axis=(1,2))
        avqv=np.mean(qv,axis=(1,2))
        press_idx=np.argwhere(avpress<qv_press*100)-1

        press700=abs(avpress-70000).argmin()
        press1000=abs(avpress-100000).argmin()
        lts=theta[press700,:,:]-theta[press1000,:,:]

        #
        # CALCULATE TCWV AND SORT 
        #
        print("calculating and sorting tcwv")
        diffhgt=np.diff(hgt,axis=0)
        rhohgt=press/(Rd*tempv) # density
        masshgt=diffhgt*rhohgt # density
        rhohgt=np.mean(rhohgt,axis=(1,2))
        tcwv=np.sum(masshgt[press_idx,:,:]*qv[press_idx,:,:],axis=0)
        print("min max tcwv ",np.min(tcwv),np.max(tcwv))
        # sort the tcwv
        isort=np.argsort(tcwv.flatten()) # sort tcwv

        # qvpert and heating rates...
        qvmean=np.mean(qv,axis=(1,2))
        qvpert=np.zeros_like(tcwv)
        lwmean=np.mean(lw_hgt,axis=(0,1))
        swmean=np.mean(sw_hgt,axis=(0,1))
        lwpert=np.zeros_like(tcwv)
        sstpert=sst.flatten()-np.mean(sst)
        lhpert=lh.flatten()-np.mean(lh)
        ltspert=lts.flatten()-np.mean(lts)

        sst_hov[irun,itime,:]=np.mean(sstpert[isort].reshape(-1,hovchunk),axis=1)
        lh_hov[irun,itime,:]=np.mean(lhpert[isort].reshape(-1,hovchunk),axis=1)
        tcwv_hov[irun,itime,:]=np.mean(tcwv.flatten()[isort].reshape(-1,hovchunk),axis=1)
        lts_hov[irun,itime,:]=np.mean(ltspert[isort].reshape(-1,hovchunk),axis=1)
        masshgt=np.mean(masshgt,axis=(1,2))

        #Pert
        for iz in range(nz):
           qvpert=100.*(qv[iz,:,:]-qvmean[iz])/qvmean[iz]
           #qratpert=(qv[iz,:,:]-qvmean[iz])/qvmean[iz]
           qlpert=ql[iz,:,:].flatten()
           qipert=qi[iz,:,:].flatten()
           qvpert=qvpert.flatten()
           lwpert=(lw_hgt[:,:]-lwmean).flatten()
           swpert=(sw_hgt[:,:]-swmean).flatten()
           wp=wvel[iz,:,:].flatten()-np.mean(wvel[iz,:,:])

           lw_hovf[irun,iz,itime,:]=np.mean(lwpert[isort].reshape(-1,hovchunk),axis=1)
           sw_hovf[irun,iz,itime,:]=np.mean(swpert[isort].reshape(-1,hovchunk),axis=1)
           qv_hov[irun,iz,itime,:]=np.mean(qvpert[isort].reshape(-1,hovchunk),axis=1)
           temp_hov[irun,iz,itime,:]=np.mean(temp[iz,:,:].flatten()[isort].reshape(-1,hovchunk),axis=1)

           ql_hov[irun,iz,itime,:]=np.mean(qlpert[isort].reshape(-1,hovchunk),axis=1)
           qi_hov[irun,iz,itime,:]=np.mean(qipert[isort].reshape(-1,hovchunk),axis=1)
           wp_hov[irun,iz,itime,:]=np.mean(wp[isort].reshape(-1,hovchunk),axis=1)
           sf_hov[irun,iz,itime,:]=np.mean(wp[isort].reshape(-1,hovchunk),axis=1)*rhohgt[iz]
           sf_hov[irun,iz,itime,:]=np.cumsum(sf_hov[irun,iz,itime,:])
        
        for iflx in flxlist:
            flxname=iflx["nm"]
            flx=np.array(ds0.variables[flxname][slice+1,:,:]).flatten()
            flx-=flx.mean() # anomaly

            # sort according to TCWV and then chunk-mean
            hov1=np.mean(flx[isort].reshape(-1,hovchunk),axis=1)

            # store to dictory slice
            #hov[flxname]=np.vstack([hov[flxname],hov1])
            hov[flxname][irun,itime,:]=hov1

       
# average across flds and plot time height
qv_hov=np.mean(qv_hov,axis=2)
ql_hov=np.mean(ql_hov,axis=2)
qi_hov=np.mean(qi_hov,axis=2)
lw_hov=np.mean(lw_hovf,axis=2)
sw_hov=np.mean(sw_hovf,axis=2)
temp_hov=np.mean(temp_hov,axis=2)
wp_hov_mean=np.mean(wp_hov,axis=2)
sf_hov_mean=np.mean(sf_hov,axis=2)

sst_av=np.mean(sst_hov,axis=1)
lh_av=np.mean(lh_hov,axis=1)
tcwv_av=np.mean(tcwv_hov,axis=1)


yticks=[1000,900,800,700,600,500,400,300,200,100]
xticks=[10,30,50,70,90]
xtickvals=[]
for irun,run in enumerate(runlist):
    idx=[bisect_left(hovx,i) for i in xticks]
    vals=[round(i,1) for i in tcwv_av[irun,idx]]
    xtickvals.append([str(i)+" ("+str(j)+")" for i,j in zip(xticks,vals)])

for irun,run in enumerate(runlist):
    print(run)
    # pick up tick vals 
    fig,(ax)=plt.subplots(nrows=2,sharex=True,figsize=(6,9))

    #
    # height plot
    #
    #pdb.set_trace()
    
    #levels=[-80,-60,-40,-20,20,40,60,80]
    levels=[-40,-30,-20,20,30,40]
    casi.make_fig(fig,ax[0],hovx,znu,qv_hov[irun,:,:],ylab="Pressure (hPa)",yticks=yticks,ymin=1000,ymax=100,cmap=cmap_r,xticks=xticks,cbar_pos=2,debug=False,conts=-999)
    csqv=ax[0].contour(hovx,znu,qv_hov[irun,:,:],colors="grey")
    plt.clabel(csqv, fontsize=10, inline=1,fmt = '%1.0f')

    ax[0].contour(hovx,znu,temp_hov[irun,:,:]-273.15,levels=[0.0],colors="purple",linewidths=2,linestyles="dashdot")

    sflevs=np.arange(-55,15,10)
    cssf=ax[0].contour(hovx,znu,100*sf_hov_mean[irun,:,:],levels=sflevs,colors="black",linewidths=2)
    plt.clabel(cssf, fontsize=10, inline=1,fmt = '%1.0f',levels=sflevs)

    qllevs=[1.e-5,3.e-5,1.e-4]
    qilevs=[1.e-6,3.e-6,1.e-5]
    csql=ax[0].contourf(hovx,znu,ql_hov[irun,:,:],levels=qllevs,hatches=["...","..","."], alpha=0,extend='max', colors = 'blue')
    csqi=ax[0].contourf(hovx,znu,qi_hov[irun,:,:],levels=qilevs,hatches=["***","**","*"], alpha=0,extend='max', colors = 'green')

    labloc="left"
    ax[0].set_title('(a)',loc=labloc)
    cs=ax[0].contour(hovx,znu,ql_hov[irun,:,:],levels=qllevs,colors="red")
    cs=ax[0].contour(hovx,znu,qi_hov[irun,:,:],levels=qilevs,colors="orange")

    ax[1].set_title('(b)',loc=labloc)
    ax[1].plot(hovx,sst_av[irun,:],color="black")
    ax[1].set_ylabel("SST (K)",color="black")
    #ax[1].set_xlabel(xlab)
    ax[1].set_xticks(list(range(0,120,20)))
    ax[1].set_ylim([-0.2,0.2])
    ax[1].tick_params(axis='y',labelcolor="black")
    ax[1].axhline(0.0,color="grey")
    ax2=ax[1].twinx()
    color='tab:red'
    ax[1].set_xlabel(xlab)
    plt.savefig(med+'/'+cvar+'/stream_'+cvar+'.jpg')
    #plt.show()
    plt.close('all')
    
    #
    #
    #

    fig,(ax)=plt.subplots(nrows=2,sharex=True)
    casi.make_fig(fig,ax[0],hovx,znu,qv_hov[irun,:,:],ylab="Pressure (hPa)",yticks=yticks,vmax=20.,ymin=1000,ymax=100,cmap=cmap_r,xticks=xticks,cbar_pos=2)
    cs=ax[0].contour(hovx,znu,wp_hov_mean[irun,:,:],[-2,-1,-0.5,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.5,1,2,5],cmap='seismic')
    ax[0].clabel(cs,inline=1,fontsize=10,fmt='%4.2f')

    ax[1].plot(hovx,sst_av[irun,:])
    ax[1].set_ylabel("SST (K)")
    ax[1].set_xlabel(xlab)
    #ax[1].set_xticks(xticks)
    #ax[1].set_xticklabels(xtickvals)

    ax2=ax[1].twinx()
    color = 'tab:red'
    ax2.set_ylabel('Latent Heat Flux (W m$^{-2}$)', color=color)  # we already handled the x-label with ax1  
    ax2.plot(hovx,lh_av[irun,:], color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.savefig(med+'/'+cvar+'/wvel_'+cvar+'.jpg')
    #plt.show()
    plt.close('all')

    #
    #
    #

    fig,(ax)=plt.subplots()
    cl=[-0.5,-0.2,0.,0.2,0.5]
    fs=9
    casi.make_fig(fig,ax,hovx,hovtimes,lts_hov[irun,:,:],title='LTS (K)',lab='',ylab=ylab,fsize=fs,cont_levs=cl,yticki=yticki)
    plt.savefig(med+'/'+cvar+'/LTS_'+cvar+'.jpg')
    #plt.show()
    plt.close('all')

    #
    #
    #
    #pdb.set_trace()    

    #Data
    ds00 = Dataset(reals+cvar+'/cloud_'+cvar+'.nc')
    #ds00 = Dataset(scra+'cloud_'+cvar+'.nc')
    cover = ds00['QICE']
    del(ds00)
    ds2 = Dataset(reals+cvar+'/HML_'+cvar+'_fldmean.nc')
    #ds2 = Dataset(scra+'HML_'+cvar+'_fldmean.nc')
    HML = ds2['HML']
    del(ds2)
    ds3 = xr.open_dataset(reals+cvar+'/u10_'+cvar+'.nc')
    ds4 = xr.open_dataset(reals+cvar+'/v10_'+cvar+'.nc')
    U10 = ds3['U10'].mean(dim=['south_north','west_east'])
    umin = ds3['U10'].min(dim=['south_north','west_east'])
    umax = ds3['U10'].max(dim=['south_north','west_east'])
    del(ds3)
    V10 = ds4['V10'].mean(dim=['south_north','west_east'])
    vmin = ds4['V10'].min(dim=['south_north','west_east'])
    vvmax = ds4['V10'].max(dim=['south_north','west_east'])
    #pdb.set_trace()
    del(ds4)
    vel = (V10*V10+U10*U10)**0.5
    ds5 = xr.open_dataset(reals+cvar+'/qv_vert_'+cvar+'.nc')
    tw = ds5.QVAPOR.mean(dim=['lat','lon'])
    twmin = ds5.QVAPOR.min(dim=['lat','lon'])
    twmax = ds5.QVAPOR.max(dim=['lat','lon'])
    sigma = ds5.QVAPOR.std(dim=['lat','lon'])
    del(ds5)
    ds6 = xr.open_dataset(reals+cvar+'/all_fldmean_'+cvar+'.nc')
    rain = ds6['RAINNC']*3600
    rain = casi.glitch(rain)
    rain = casi.glitch(rain)
    rain = rain.where(rain>0)
    LH = ds6['ACLHF']
    LH = casi.glitch(LH)
    LH = casi.glitch(LH)
    LH = LH.where(LH>50)
    SH = ds6['ACHFX']
    SH = casi.glitch(SH)
    SH = casi.glitch(SH)
    SH = SH.where(SH>4)
    LW = ds6['ACLWTNET'] - ds6['ACLWTNETC']
    LW = LW.where(LW>0)
    LW = casi.glitch(LW)
    LWC = ds6['ACLWTNETC']
    LWC = LWC.where(LWC<0)
    LWC = casi.glitch(LWC)
    SW = ds6['ACSWTNET'] - ds6['ACSWTNETC']
    SW = casi.glitch(SW)
    SW = SW.rolling(XTIME=24).mean()
    #SW = SW.where(SW>-110) 
    SWC = ds6['ACSWTNETC']
    SWC = casi.glitch(SWC)
    SWC = SWC[:,0,0].rolling(XTIME=24).mean()
    SWC = SWC.where(SWC>380)
    SWC = SWC.where(SWC<393)
    del(ds6)
    ds7 = xr.open_dataset(reals+cvar+'/TSK_'+cvar+'.nc')
    tskme = ds7.TSK.mean(dim=['south_north','west_east'])
    tskmin = ds7.TSK.min(dim=['south_north','west_east'])
    tskmax = ds7.TSK.max(dim=['south_north','west_east'])
    del(ds7)
    #pdb.set_trace()

    ##Plot
    cl=[0.0]
    fs=9
    vmax = 0.25
    limi = 55
    fig,(ax)=plt.subplots(nrows=2,ncols=5,figsize=(16,18)) #width,height
    #fig,(ax)=plt.subplots()
    plt.subplots_adjust(hspace=0.20,wspace=0.25,top=0.9,right=0.98,left=0.05)
    casi.make_fig(fig,ax[0,0],hovx,hovtimes,sst_hov[irun,:,:],lab='',ylab=ylab,fsize=fs,cont_levs=-999,yticki=yticki,vmax=vmax)
    ax[0,0].set_xlabel('TCWV %-tile',fontsize=fs)
    ### Cloud
    ax[0,1].plot(cover[0:len(hovtimes)+0,0,0]*100,hovtimes,color='k')
    ax[0,1].set_ylim(0,limi)
    if cvar == 'uu':
        ax[0,1].set_xlim(25,100)
    ax[0,1].set_xlabel('Cloud cover (%)',fontsize=fs)
    #ax[0,1].set_yticks([0,4,8,16,20,24])
    ax[0,1].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ## HML
    ax[0,2].plot(HML[0:len(hovtimes),0,0],hovtimes,color='darkblue')
    ax[0,2].set_ylim(0,limi)
    if cvar == 'uu':
        ax[0,2].set_xlim(3,12) #3,10
    ax[0,2].set_xlabel('Mixed Layer Depth (m)',fontsize=fs)
    #ax[0,2].set_yticks([0,4,8,16,20,24])
    ax[0,2].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ### Wind
    ax[0,3].plot(U10[0:len(hovtimes)],hovtimes,color='b',label='U')
    ax[0,3].plot(V10[0:len(hovtimes)],hovtimes,color='r',label='V')
    ax[0,3].plot(vel[0:len(hovtimes)],hovtimes,color='purple',label='Speed')
    ax[0,3].set_ylim(0,limi)
    #ax[0,3].set_xlim(-10,10)
    ax[0,3].set_xlabel('10m - Wind (m/s)',fontsize=fs)
    #ax[0,3].set_yticks([0,4,8,16,20,24])
    ax[0,3].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ax[0,3].legend(loc='upper center',bbox_to_anchor=(0.5, 1.10),fontsize=fs,ncol=3)
    ###TCWV
    ax[0,4].plot(twmax[0:len(hovtimes)],hovtimes,color='b',label = 'Max')
    ax[0,4].plot(tw[0:len(hovtimes)],hovtimes,color='green',label = 'Mean')
    ax[0,4].plot(twmin[0:len(hovtimes)],hovtimes,color='r',label = 'Min')    
    ax[0,4].set_ylim(0,limi)
    #ax[0,4].set_xlim(15,90)
    ax[0,4].set_xlabel('TPW (mm)',fontsize=fs)
    #ax[0,4].set_yticks([0,4,8,16,20,24])
    ax[0,4].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ax[0,4].legend(loc='upper center',bbox_to_anchor=(0.5, 1.1),fontsize=fs,ncol=3)
    ### Hov muller 
    casi.make_fig(fig,ax[1,0],hovx,hovtimes,sst_hov[irun,:,:],lab='',ylab=ylab,fsize=fs,cont_levs=-999,yticki=yticki,vmax=vmax)
    ax[1,0].set_xlabel('TCWV %-tile',fontsize=fs)
    #sigma
    ax[1,1].plot(sigma[0:len(hovtimes)],hovtimes, color = 'royalblue')
    ax[1,1].set_ylim(0,limi)
    if cvar == 'qv_time':
         ax[1,1].set_xlim(2,5.5) 
    #ax[1,1].set_xlim(2,16)
    ax[1,1].set_xlabel('$\sigma_{TPW}$ (mm)',fontsize=fs)
    #ax[1,1].set_yticks([0,4,8,16,20,24])
    ax[1,1].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ##Rain     
    ax[1,2].plot(rain[0:len(hovtimes),0,0],hovtimes,color='blue')
    ax[1,2].set_ylim(0,limi)
    if cvar == 'qv_time':
        ax[1,2].set_xlim(0,0.6) #0,1
    ax[1,2].set_xlabel('Mean Precip (mm/hr)',fontsize=fs)
    #ax[1,2].set_yticks([0,4,8,16,20,24])
    ax[1,2].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #pdb.set_trace()
    ##Surface fluxes
    #Latent heat
    ax[1,3].plot(LH[0:len(hovtimes),0,0],hovtimes,color='darkcyan')
    ax[1,3].set_ylim(0,limi)
    #ax[1,3].set_xlim(80,220) #80,220
    ax[1,3].set_xlabel('Latent heat (W m$^{-2}$)',fontsize=fs)
    #ax[1,3].set_yticks([0,4,8,16,20,24])
    ax[1,3].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #Sensible heat 
    ax[1,4].plot(SH[0:len(hovtimes),0,0],hovtimes,color='darkcyan')
    ax[1,4].set_ylim(0,limi)
    if cvar == 'qv_time':
        ax[1,4].set_xlim(3,12) #3,15
    ax[1,4].set_xlabel('Sensible heat (W m$^{-2}$)',fontsize=fs)
    #ax[1,4].set_yticks([0,4,8,16,20,24])
    ax[1,4].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    plt.savefig(med+'/'+cvar+'/sst_hov_'+cvar+'_one.jpg', bbox_inches = 'tight')
    plt.savefig(med+'/'+cvar+'/sst_hov_'+cvar+'_one.pdf', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')
    #############
    fig,(ax)=plt.subplots(nrows=2,ncols=5,figsize=(16,18)) #width,height
    #fig,(ax)=plt.subplots()
    plt.subplots_adjust(hspace=0.25,wspace=0.25,top=0.9,right=0.98,left=0.05)
    ### Radiation
    ### Hov muller
    casi.make_fig(fig,ax[0,0],hovx,hovtimes,sst_hov[irun,:,:],lab='',ylab=ylab,fsize=fs,cont_levs=-999,yticki=yticki,vmax=vmax)
    ax[0,0].set_xlabel('TCWV %-tile',fontsize=fs)
    #LW
    ax[0,1].plot(LW[0:len(hovtimes),0,0],hovtimes,color='r')
    ax[0,1].set_ylim(0,limi)
    #ax[0,1].set_xlim(15,125) #15,125
    ax[0,1].set_xlabel('LW CRE (W m$^{-2}$)',fontsize=fs)
    #ax[0,1].set_yticks([0,4,8,16,20,24])
    ax[0,1].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #SW
    ax[0,2].plot(SW[0:len(hovtimes),0,0].rolling(XTIME=24).mean(),hovtimes,color='darkorange')
    ax[0,2].set_ylim(0,limi)
    #ax[0,2].set_xlim(-165,-30) #-140,-30
    ax[0,2].set_xlabel('SW CRE (W m$^{-2}$)',fontsize=fs)
    #ax[0,2].set_yticks([0,4,8,16,20,24])
    ax[0,2].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #LW C  
    ax[0,3].plot(LWC[0:len(hovtimes),0,0],hovtimes,color='r')
    ax[0,3].set_ylim(0,limi)
    #ax[0,3].set_xlim(-300,-265) #-300, -265
    ax[0,3].set_xlabel('LW Clear-Sky (W m$^{-2}$)',fontsize=fs)
    #ax[0,3].set_yticks([0,4,8,16,20,24])
    ax[0,3].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #SWC
    #pdb.set_trace()
    ax[0,4].plot(SWC[0:len(hovtimes)],hovtimes,color='darkorange')
    ax[0,4].set_ylim(0,limi)
    #ax[0,4].set_xlim(370,395)
    ax[0,4].set_xlabel('SW Clear-Sky (W m$^{-2}$)',fontsize=fs)
    #ax[0,4].set_yticks([0,4,8,16,20,24])
    ax[0,4].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ### Hov muller
    casi.make_fig(fig,ax[1,0],hovx,hovtimes,sst_hov[irun,:,:],lab='',ylab=ylab,fsize=fs,cont_levs=-999,yticki=yticki,vmax=vmax)
    ax[1,0].set_xlabel('TCWV %-tile',fontsize=fs)
    #NET TOA CRE
    ax[1,1].plot(SW[0:len(hovtimes),0,0].rolling(XTIME=24).mean()+LW[0:len(hovtimes),0,0],hovtimes,color='r')
    ax[1,1].set_ylim(0,limi)
    #ax[1,1].set_xlim(-60,50) #15,125
    ax[1,1].set_xlabel('NET CRE (W m$^{-2}$)',fontsize=fs)
    #ax[1,1].set_yticks([0,4,8,16,20,24])
    ax[1,1].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #NET TOA Clear
    ax[1,2].plot(SWC[0:len(hovtimes)]+LWC[0:len(hovtimes),0,0],hovtimes,color='r')
    ax[1,2].set_ylim(0,limi)
    #ax[1,2].set_xlim(90,130) 
    ax[1,2].set_xlabel('NET Clear (W m$^{-2}$)',fontsize=fs)
    #ax[1,2].set_yticks([0,4,8,16,20,24])
    ax[1,2].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    #U
    ax[1,3].plot(umin[0:len(hovtimes)],hovtimes,color='r',label='Min')
    ax[1,3].plot(U10[0:len(hovtimes)],hovtimes,color='green',label='Mean')
    ax[1,3].plot(umax[0:len(hovtimes)],hovtimes,color='b',label='Max')
    ax[1,3].set_ylim(0,limi)
    #ax[1,3].set_xlim(-25,20)
    ax[1,3].set_xlabel('u-10m (m/s)',fontsize=fs)
    #ax[1,3].set_yticks([0,4,8,16,20,24])
    ax[1,3].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    ax[1,3].legend(loc='upper center',bbox_to_anchor=(0.5, 1.1),fontsize=fs,ncol=3)
    #V
    #pdb.set_trace()
    ax[1,4].plot(vmin[0:len(hovtimes)],hovtimes,color='r',label='v-min')
    ax[1,4].plot(V10[0:len(hovtimes)],hovtimes,color='green',label='v-mean')
    ax[1,4].plot(vvmax[0:len(hovtimes)],hovtimes,color='b',label='v-max')
    ax[1,4].set_ylim(0,limi)
    #ax[1,4].set_xlim(-20,20)
    ax[1,4].set_xlabel('v-10m (m/s)',fontsize=fs)
    #ax[1,4].set_yticks([0,4,8,16,20,24])
    ax[1,4].set_yticks([0,5,10,15,20,25,30,35,40,45,50,55])
    plt.savefig(med+'/'+cvar+'/sst_hov_'+cvar+'_two.jpg', bbox_inches = 'tight')
    plt.savefig(med+'/'+cvar+'/sst_hov_'+cvar+'_two.pdf', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')


hov["ATM_TOTAL"]=np.zeros([nrun,ntime,nhov])
hov["SFC_TOTAL"]=np.zeros([nrun,ntime,nhov])
hov["TOA_TOTAL"]=np.zeros([nrun,ntime,nhov])
hov["TOA_CLEAR"]=np.zeros([nrun,ntime,nhov])
hov["TOA_CRE"]=np.zeros([nrun,ntime,nhov])
    
for irun,run in enumerate(runlist):
    print("almagated plots",run)
    # Hov is done, so just need to plot, oo yes.
    # SPECIAL PLOT incl CLEAR and CLD 

    fig,(ax)=plt.subplots(ncols=3,nrows=2,figsize=(6,6),sharey=True)
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.2,wspace=0.2)

    fs=9

    yticki=5
    cl=[-30,-20,-10,10,20,30]

    vmax=30.
    xticks=[0,25,50,75,100]

    flxname="ACLHF"
    iflx=next(flx for flx in flxlist if flx["nm"]==flxname)
    casi.make_fig(fig,ax[0,0],hovx,hovtimes,hov[flxname][irun,:,:],title=iflx['title'],lab='(a)',ylab=ylab,fsize=fs,vmax=vmax,cbar_pos=0,xticks=xticks,yticki=yticki)

    flxname="ACHFX"
    iflx=next(flx for flx in flxlist if flx["nm"]==flxname)
    casi.make_fig(fig,ax[1,0],hovx,hovtimes,hov[flxname][irun,:,:],title=iflx['title'],lab='(b)',xlab=xlab,ylab=ylab,fsize=fs,vmax=vmax,cbar_pos=0,xticks=xticks,yticki=yticki)

    # clear LW
    flxname="ACLWTBNETC"
    casi.make_fig(fig,ax[0,1],hovx,hovtimes,hov[flxname][irun,:,:],title="LW Clear",lab='(c)',kfac=1,fsize=fs,vmax=vmax,cbar_pos=0,xticks=xticks,yticki=yticki)
    # cloudy
    flxname="ACLWTBNET"
    casi.make_fig(fig,ax[1,1],hovx,hovtimes,hov["ACLWTBNET"][irun,:,:]-hov["ACLWTBNETC"][irun,:,:],title="LW Cloud",lab='(d)',kfac=1,fsize=fs,xlab=xlab,vmax=vmax,cbar_pos=0,xticks=xticks,yticki=yticki)

    flxname="ACSWTBNETC"
    casi.make_fig(fig,ax[0,2],hovx,hovtimes,hov[flxname][irun,:,:],title='SW Clear',lab='(e)',fsize=fs,vmax=vmax,xticks=xticks,yticki=yticki)

    flxname="ACSWTBNET"
    casi.make_fig(fig,ax[1,2],hovx,hovtimes,hov["ACSWTBNET"][irun,:,:]-hov["ACSWTBNETC"][irun,:,:],title='SW Cloud ',lab='(f)',xlab=xlab,fsize=fs,vmax=vmax,xticks=xticks,yticki=yticki)
    plt.savefig(med+'/'+cvar+'/hovmoller_fluxes_'+cvar+'.jpg')
    #plt.show()
    plt.close('all')
    
    # total atmosphere and total surface append:
    flxlist.append({"nm":"ATM_TOTAL","title":"total atmospheric heating anomaly"})
    hov["ATM_TOTAL"][irun,:,:]=hov["ACLHF"][irun,:,:]+hov["ACHFX"][irun,:,:]+hov["ACLWTBNET"][irun,:,:]+hov["ACSWTBNET"][irun,:,:]
    flxlist.append({"nm":"SFC_TOTAL","title":"surface anomaly"})
    hov["SFC_TOTAL"][irun,:,:]=hov["ACLWBNET"][irun,:,:]+hov["ACSWBNET"][irun,:,:]+hov["ACLHF"][irun,:,:]+hov["ACHFX"][irun,:,:]
    hov["TOA_TOTAL"][irun,:,:]=hov["ACLWTNET"][irun,:,:]+hov["ACSWTNET"][irun,:,:]
    hov["TOA_CLEAR"][irun,:,:]=hov["ACLWTNETC"][irun,:,:]+hov["ACSWTNETC"][irun,:,:] 
    hov["TOA_CRE"][irun,:,:]=hov["TOA_TOTAL"][irun,:,:]+hov["TOA_CLEAR"][irun,:,:]
    
    vmax=50
    # total plots 
    fig,(ax)=plt.subplots(figsize=(10,6),ncols=2)
    casi.make_fig(fig,ax[0],hovx,hovtimes,hov["ATM_TOTAL"][irun,:,:],title="total atmos",lab="(a)",xlab=xlab,ylab=ylab,yticki=yticki,vmax=vmax)   
    casi.make_fig(fig,ax[1],hovx,hovtimes,hov["SFC_TOTAL"][irun,:,:],title="total surface",lab="(b)",xlab=xlab,yticki=yticki,vmax=vmax) 
    plt.savefig(med+'/'+cvar+'/total_'+cvar+'.jpg')
    #plt.show()
    plt.close('all')
    
    vmax=80
    fig,(ax)=plt.subplots(figsize=(15,6),ncols=3)
    casi.make_fig(fig,ax[0],hovx,hovtimes,hov["TOA_TOTAL"][irun,:,:],title="NTOA all-sky",lab="(a)",xlab=xlab,ylab=ylab,yticki=yticki,vmax=vmax)
    casi.make_fig(fig,ax[1],hovx,hovtimes,hov["TOA_CRE"][irun,:,:],title="NTOA CRE",lab="(b)",xlab=xlab,yticki=yticki,vmax=vmax)
    casi.make_fig(fig,ax[2],hovx,hovtimes,hov["TOA_CLEAR"][irun,:,:],title="NTOA Clear",lab="(c)",xlab=xlab,yticki=yticki,vmax=40)
    plt.savefig(med+'/'+cvar+'/TOA_'+cvar+'.jpg')
    #plt.show()
    plt.close('all')

    #pdb.set_trace()

print('#####Complete!#####')

