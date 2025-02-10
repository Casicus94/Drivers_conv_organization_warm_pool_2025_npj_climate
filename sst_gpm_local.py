
import pandas as pd
from netCDF4 import Dataset, num2date, date2index
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as collections
import matplotlib.colorbar as cbar
import matplotlib.gridspec as gridspec
import numpy as np
#import mdir 
import sys
import timeit
import xarray as xr
import os
from scipy import stats
import dateutil.parser
from datetime import datetime
from getargs import getargs
import setup 
from string import ascii_lowercase

from dir_locs import locs
import pdb

# pdb.set_trace() # WOW!

#
# still to do 
#
# 1. replace min, mean and max of W with block plot and add lines to that for mean etc. 
# 2. delete sigma W
# 3. insert plot of mixed layer depth derived from Hinawari SST diurnal cycle.
# 4. insert plot of OLR *IF* we can find a reasonable product


#dir="/afs/ictp.it/home/m/mde_vera/temporary/"
#dir="/media/mde_vera/Seagate Backup Plus Drive/w_5_15_data/"
#dir=r"D:\w_5_15_data/"
dir=r"D:\2_9_data/"

utc=9

#composite settings:
nextrad=3 ### NUMBER OF DAYS TO SAMPLE PRIOR/AFTER THE REVERSAL 
# size of interpolated events is constant... we can save in a dictionary of 2d numpy arrays to average at the end ?
ninterpol=201 ### NUMBER OF INTERPOLATION POINTS USED

def add_rev(ax,fmin,fmax,reverses):
    for rev in reverses:
        #ax[2].axhline(y=rev["start"], linewidth=lwd, linestyle=':', color='red')
        #ax[2].axhline(y=rev["end"], linewidth=lwd, linestyle=':', color='green')
        ax.add_patch(patches.Rectangle(
               xy=(fmin,rev["start"]),  # point of origin.
               width=fmax-fmin, height=rev["end"]-rev["start"], linewidth=1,
               color='red', fill=True, alpha=0.15))

def add_zero(hov,hov1):
    for key in hov.keys():
        hov[key]=np.vstack([hov[key],hov1[key]])
    return(hov)

def oneplot(date1,pars):
    # composite settings:
    nextra=pars["nextrad"]*24 # number of hourly steps to grab before start and after end of rev event
    plt_type=".pdf"

    region="wpac_2-9N"
    date2=date1+pd.DateOffset(days=pars["dateinc"], hours=23)
    daterangeday = pd.date_range(date1,date2)
    #print(daterange)
    
    #subtract UTC time to get the timestamp of data
    #locdate1=date1+pd.DateOffset(hours=-utc)
    #locdate2=date2+pd.DateOffset(hours=-utc)
    daterangehour = pd.date_range(date1,date2, freq='H')
    #print(daterange)

    ntime=0 # zero row at start

    # set up SST bins and sum vector

    tpwb={'min':25,'max':85,'inc':5}
    sstb={'min':25,'max':33,'inc':0.2}
    sstamin=-1;sstamax=1
    sstmin=27;sstmax=31

    npts=70*100 # hardwired for moment
    hovchunk=70
    nhov=int(npts/hovchunk)
    hovx=100*(np.arange(nhov)+0.5)/nhov

    bnds={'tpw':tpwb,'sst':sstb}
    vartitle={'tpw':"tcwv (kg/m2)",'sst':'SST (oC)','gpm':'precipitation (mm day$^{-1}$)'}

    vars=bnds.keys()
    bins,nbins,ntot,ptot,npostot,ppostot,p2tot,hov,hovmax,hovmin={},{},{},{},{},{},{},{},{},{}

    for var in vars:
        n=int((bnds[var]["max"]-bnds[var]["min"])/bnds[var]["inc"]+1)
        #print ( bnds[var]["min"],bnds[var]["max"],var)
        bins[var]=np.linspace(int(bnds[var]["min"]),int(bnds[var]["max"]),n)
        nbins[var]=len(bins[var])-1
        ntot[var]=np.zeros(nbins[var],dtype=np.int64)
        ptot[var]=np.zeros(nbins[var],dtype=np.longdouble)
        npostot[var]=np.zeros(nbins[var],dtype=np.int64)
        ppostot[var]=np.zeros(nbins[var],dtype=np.longdouble)
        p2tot[var]=np.zeros(nbins[var],dtype=np.longdouble)
        hov[var]=np.zeros(nhov)
        hovmax[var]=np.zeros(nhov)
        hovmin[var]=np.zeros(nhov)

    hov['precip']=np.zeros(nhov)

 

    ytick=[]
    yticklabels=[]
    ytickcycle=int(len(daterangeday)/10) # 10 labels...
    #print(len(locdaterange))

#    meanflds={[],[],[]}
    # this should be in a dictionary...


    # this is really yukky, replace by pandas dataframe
    minsst=[]
    maxsst=[]
    meansst=[]
    stdevW=[]
    minW=[]
    maxW=[]
    meanW=[]
    sst_hires_mean=[]
    sst_hires_max=[]
    sst_hires_min=[]
    mld=[]
    diff=[]
    slopes=[]
    slopes_np=[]
    years=[]
    months=[]
    days=[]
    hours=[]
    olr_clr=[]
    olr_all=[]
    sw_clr=[]
    sw_all=[]
    minu=[]
    maxu=[]
    avgu=[]
    mjostrength=[]
    mjophase=[]
    meanpp=[]

    df=pd.DataFrame()

    print (str(date1))
    time_units="hours since "+str(date1)

    ####
    # SLOPES FILES
    ####
    fileout = Dataset('./slopes.nc', 'w')
    fileout.createDimension('time', None)
    times_out = fileout.createVariable('time', np.float32, 'time',)
    times_out.setncattr('units',time_units)
    times_out[:] = np.arange(len(daterangehour))


    for idate,sdate in enumerate(daterangehour):
        fulldate=sdate.strftime("%Y%m%d%H")
        yy=sdate.strftime("%Y")
        mm=sdate.strftime("%m")
        dd=sdate.strftime("%d")
        hh=sdate.strftime("%H")
        
        date=yy+mm+dd
        tag="_"+region+"_"+fulldate+".nc"
        tagdaily="_"+region+"_"+date+".nc"
        tagdailymean="_"+region+"_"+date+"_mean.nc"

        print(sdate)
        # OLR
        if False:
            try:
                dolr=Dataset(sdir["wrk"]+'olr'+tagdailymean) 
                rad=dolr.variables["toa_lw_clr_1h"][hh,0,0]
                olr_clr.append(rad)
                rad=dolr.variables["toa_lw_all_1h"][hh,0,0]
                olr_all.append(rad)
                rad=dolr.variables["toa_sw_clr_1h"][hh,0,0]
                sw_clr.append(rad)
                rad=dolr.variables["toa_sw_all_1h"][hh,0,0]
                sw_all.append(rad)
                dolr.close()
            except:
                print ("missing olr file",sdir["wrk"]+'olr'+tag)
                olr_clr.append(np.nan)
                olr_all.append(np.nan)
                sw_clr.append(np.nan)
                sw_all.append(np.nan)

        # ERA5 RADN - ALREADY THE MEAN
        try:
            erafile=sdir["wrk"]+"tropbox_radn_"+yy+"_mean.nc"
            dolr=xr.open_dataset(erafile)
            tim=datetime(int(yy),int(mm),int(dd),int(hh))
            fluxes=dolr.sel(time=slice(tim,tim))
            olr_all.append(np.squeeze(np.array(fluxes.variables["ttr"])))
            sw_all.append(np.squeeze(np.array(fluxes.variables["tsr"])))
            olr_clr.append(np.squeeze(np.array(fluxes.variables["ttrc"])))
            sw_clr.append(np.squeeze(np.array(fluxes.variables["tsrc"])))
            dolr.close()
        except:
            print("something went wrong with ERA radn")
            exit()
 

        # tick per day with date
        if np.mod(idate,ytickcycle)==0:
            localdate=sdate+pd.DateOffset(hours=utc)
            #print(localdate.strftime("%d%b%Y"))
            ytick.append(ntime)
            yticklabels.append(localdate.strftime("%d%b%Y"))

        #sstfile=mdir.scr+"obs/sst_"+region+"_"+date+".nc"
        sstfile=sdir["wrk"]+'sst_'+region+"_"+date+".nc"
        try:
            dsst=Dataset(sstfile) 
        except:
            print("file error here:",sstfile)
        
        sst=dsst.variables["sst"]
        lvar={}
        lvar['sst']=np.array(sst).flatten()
        dsst.close()

        # place row of zeros in hov, to be overwritten if file there
        hov1,hov1max,hov1min={},{},{}
        for key in hov.keys():
            hov1[key]=np.zeros(nhov)
            hov1max[key]=np.zeros(nhov)
            hov1min[key]=np.zeros(nhov)

        # himawari SST
        sstfile2=sdir["wrk"]+'sst_hires'+tag
        try:
            dsst=Dataset(sstfile2)
            sst_hires=dsst.variables["sea_surface_temperature"][:]
            idx=np.where(sst_hires.mask==True)
            sst_hires[idx]=np.nan
# nan check here
            try:
              sst_hires_min.append(np.nanmin(sst_hires))  ##### This has an error,
              sst_hires_max.append(np.nanmax(sst_hires))
              sst_hires_mean.append(np.nanmean(sst_hires))
            except:
              print("file error:",sstfile2)
              continue
              
            dsst.close()
        except:
            print("file missing:",sstfile2)
            sst_hires[:]=np.nan
            sst_hires_min.append(np.nan)
            sst_hires_max.append(np.nan)
            sst_hires_mean.append(np.nan)
       
        # mimic TCWV
        try:
            #dmim=Dataset(mdir.scr+"obs/mim"+tag) 
            dmim=Dataset(sdir["wrk"]+'mim'+tag) 
            tpw=dmim.variables["tpwGrid"]
            minW.append(np.nanmin(tpw))
            maxW.append(np.nanmax(tpw))
            meanW.append(np.nanmean(tpw))
        except:
            #print("error opening file",mdir.scr+"obs/mim"+tag)
            print("error opening file",sdir["wrk"]+'mim'+tag)
            minW.append(np.nan)
            maxW.append(np.nan)
            meanW.append(np.nan)
            #hov=add_zero(hov,hov1)
            #hovmax=add_zero(hovmax,hov1)
            #hovmin=add_zero(hovmin,hov1)
        #print(ntime+1,len(minW))
        if ntime+1!=len(minW):
            exit("error len")

        #if "tpwGrid" in dmim.variables:
        tpw=dmim.variables["tpwGrid"]
        lvar['tpw']=np.array(tpw).flatten() # sst or tcwv
        dmim.close()

        # Himawari based MLD - THIS IS DAILY...
        mldfile=sdir["wrk"]+'mld'+tagdaily
        try:
            dmld=Dataset(mldfile) 
            mld.append(np.array(dmld.variables["mld"]).flatten())
            dmld.close()
        except:
            #print("error opening file",mdir.scr+"obs/mim"+tag)
            print("error opening file: ",mldfile)
            mld.append(np.nan)
        
        #else:
        #    lvar['tpw']=np.zeros(10000)
            
        # GPM precip
        try:
            #dgpm=Dataset(mdir.scr+"obs/gpm"+tag) 
            dgpm=Dataset(sdir["wrk"]+'gpm'+tag)
            precip=dgpm.variables["precipitationCal"]
            precip=np.array(precip).flatten()
            meanpp.append(np.nanmean(precip)) # SHOULD USE COSINE WEIGHTING.
            dgpm.close()
        except:
            #print("error opening file:",mdir.scr+"obs/gpm"+tag)
            print("error opening file:",sdir["wrk"]+'gpm'+tag)
            #hov=add_zero(hov,hov1)
            #hovmax=add_zero(hovmax,hov1)
            #hovmin=add_zero(hovmin,hov1)
            meanpp.append(np.nan)

        #standard deviation of W
        try:
            #dmim=Dataset(mdir.scr+"obs/mim"+tag) 
            dstdW=Dataset(sdir["wrk"]+'stdW'+tag)
            stdevW.append(dstdW.variables["tpwGrid"][:][0]) 
        except:
            #print("error opening file",mdir.scr+"obs/mim"+tag)
            print("error opening file",sdir["wrk"]+'stdW'+tag)
            #stdevW.append(np.nan)

        avguwindfile=sdir["wrk"]+'fldmean_10m_u_'+region+"_"+yy+".nc"
        try:
            duwind=Dataset(avguwindfile) 
        except:
            print("file error:",avguwindfile)
        
        all_times=duwind.variables["time"]
        timenow=sdate.strftime("%Y-%m-%dT%H:00:00")
        #print(timenow)

        dt = dateutil.parser.parse(timenow)
        dt_idx = date2index(dt, all_times)        
        avgu.append(duwind.variables['u10'][dt_idx,0,0])
        uwindfile=sdir["wrk"]+'10m_u_'+region+"_"+yy+".nc"
        try:
            uwind=Dataset(uwindfile) 
        except:
            print("file error:",uwindfile)
        
        minu.append(np.nanmin(uwind.variables['u10'][dt_idx]))
        maxu.append(np.nanmax(uwind.variables['u10'][dt_idx]))

        #mjo
        omi = pd.read_csv(sdir["scripts"]+"omi.csv")
        omi.columns= ["year","month","day","hour","PC1","PC2","PC1+2"]
        
        mjo_amp=omi.loc[(omi['year'] == int(yy)) & (omi['month'] == int(mm)) & (omi['day'] == int(dd)), "PC1+2"].item()
        mjostrength.append(mjo_amp)
        
        mjox=omi.loc[(omi['year'] == int(yy)) & (omi['month'] == int(mm)) & (omi['day'] == int(dd)), "PC2"].item()
        mjoy=-omi.loc[(omi['year'] == int(yy)) & (omi['month'] == int(mm)) & (omi['day'] == int(dd)), "PC1"].item()
        
        if((mjox > 0) & (mjoy > 0)):
            if(mjox < mjoy):
                mjophase.append(6)
            else:
                mjophase.append(5)
        elif((mjox < 0) & (mjoy > 0)):
            if(abs(mjox) < mjoy):
                mjophase.append(7)
            else:
                mjophase.append(8)
        elif((mjox < 0) & (mjoy < 0)):
            if(mjox < mjoy):
                mjophase.append(1)
            else:
                mjophase.append(2)
        else:
            if(mjox < abs(mjoy)):
                mjophase.append(3)
            else:
                mjophase.append(4)
        
        #min max SST
        minsst.append(np.nanmin(lvar['sst']))
        maxsst.append(np.nanmax(lvar['sst']))
        
        # hov stuff
        isort=np.argsort(lvar['tpw']) # sort tcwv
        
        # this is the clever part where we sort TCWV and then chunk all data by hovchunk.
        for var in vars:
            hov1[var]=np.mean(np.array(lvar[var])[isort].reshape(-1,hovchunk),axis=1)
            if var=='sst':
                msst=np.mean(hov1[var])
                meansst.append(msst)
                hov1[var]-=msst # SST anomaly about slice mean...
            hov1max[var]=np.max(np.array(lvar[var])[isort].reshape(-1,hovchunk),axis=1)
            hov1min[var]=np.min(np.array(lvar[var])[isort].reshape(-1,hovchunk),axis=1)
            #print(hov[var].shape)
            #print(hov1[var].shape)
            hov[var]=np.vstack([hov[var],hov1[var]])
            hovmax[var]=np.vstack([hovmax[var],hov1max[var]])
            hovmin[var]=np.vstack([hovmin[var],hov1min[var]])
        # yuk
        hov['precip']=np.vstack([hov['precip'],np.mean(np.array(precip)[isort].reshape(-1,hovchunk),axis=1)])

        ntime+=1 
        years.append(yy)
        months.append(mm)
        days.append(dd)
        hours.append(hh)

        idx=np.argwhere(precip>0.01)

        for var in vars:
            n, _   = np.histogram(lvar[var], bins=bins[var])
            sy, _  = np.histogram(lvar[var], bins=bins[var], weights=precip)
            sy2, _ = np.histogram(lvar[var], bins=bins[var], weights=precip*precip)
            npos, _  = np.histogram(lvar[var][idx], bins=bins[var])
            sypos, _  = np.histogram(lvar[var][idx], bins=bins[var], weights=precip[idx])
            ntot[var]+=n
            ptot[var]+=sy
            npostot[var]+=npos
            ppostot[var]+=sypos
            p2tot[var]+=sy2


    #----------------------END OF LONG TIME LOOP

    # flip the OLR 
    olr_all=[-x for x in olr_all]
    olr_clr=[-x for x in olr_clr]


    #   
    # calculate slopes
    #
    x = np.array(hovx)
    for i in range(1, ntime+1):
        y = np.array(hov['sst'][i])
        xm = np.ma.masked_array(x,mask=np.isnan(y)).compressed()
        ym = np.ma.masked_array(y,mask=np.isnan(y)).compressed()

        ynp=y.copy()
        ynp=np.where(hov['precip'][i]<1.0,ynp,np.nan)
        
        xnpm=np.ma.masked_array(x,mask=np.isnan(ynp)).compressed()
        ynpm=np.ma.masked_array(ynp,mask=np.isnan(ynp)).compressed()

        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(xm, ym)
            slopes.append(slope)  # multiply by 100 to give K as units (K %-1)
        except:
            print("Error: ",i)
            slopes.append(np.nan)

        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(xnpm, ynpm)
            slopes_np.append(slope)  # multiply by 100 to give K as units (K %-1)
        except:
            print("Error np: ",i)
            slopes_np.append(np.nan)

    #pdb.set_trace()

    # negative standard deviation 
    slopes_smth=pd.Series(slopes).rolling(window=24,min_periods=1,center=True).mean()
    slopes_np_smth=pd.Series(slopes_np).rolling(window=24,min_periods=1,center=True).mean()

    slopes_av=np.mean(slopes_smth)
    slopes_sd=np.std(slopes_smth)
    fac=0.3
    slopes_crit_lo=-fac*slopes_sd # slopes_sd # critical threshold for an event
    slopes_crit_hi=fac*slopes_sd # slopes_sd # critical threshold for an event
     
    #------------------------------
    # catalogue all reversal events
    #------------------------------
    times=np.arange(ntime)

    # make series of 0=normal, 1=reversal #### NOTE CRITICAL THRESHOLD HERE < -500 switches to use sigma

    ###
    ### NEW CONDITION - SLOPES ARE NEGATIVE AND THE MINMUM < 2*sd

    slopes_crit_lo=pars["event_thresh"] if pars["event_thresh"]>-500 else slopes_crit_lo
    slopes_crit_hi=-pars["event_thresh"] if pars["event_thresh"]>-500 else slopes_crit_hi
    # slopes_mask=np.where(np.array(slopes_smth)<slopes_crit_lo,1,0)
    slopes_mask=np.where(np.array(slopes_smth)<0.0,1,0)
    slopes_mask_hi=np.where(np.array(slopes_smth)>slopes_crit_hi,1,0)

    # running sum, will equal 1 when switching from normal from/to reversal
    slopes_idx_run=pd.Series(slopes_mask).rolling(window=2,min_periods=1).sum().astype(int)

    # pick up switches between states
    slopes_idx=np.where(slopes_idx_run==1)[0]

    # safety: if we started already in a reversal, throw first case away
    if slopes_mask[0]==1:
        print ("REJECTING FIRST POINT",len(slopes_idx))
        print (slopes_idx_run[0:10])
        print ("smth ",slopes_smth[0:10])
        slopes_idx=slopes_idx[2:] # THIS WAS A BUG - need to throw both first entries away!

    # safety: if we end in a reversal (list is odd length), throw last case away 
    if len(slopes_idx)%2==1:
        slopes_idx=slopes_idx[:-1]

    # construct a list of reversals
    reverses=[]

    # start and ends are even and odd entries:
    r_starts=slopes_idx[::2]
    r_ends=slopes_idx[1::2]
    r_lens=r_ends-r_starts+1         # length of events 
    r_gaps=r_starts[1:]-r_ends[:-1]  # gap between events (#1 less than lens)
    print("starts:",r_starts)
    print("ends:",r_ends)
    print("lens:",r_lens)
    print("gaps:",r_gaps)
    # join events together if the gap is smaller than the minimum length
    gdel=[]
    for ig,g in enumerate(r_gaps):
        if g<pars["gap_len"]:
            gdel.append(ig)
            print (ig, " too short") 

    # join up series here
    if len(gdel)>0:
        r_ends=np.delete(r_ends,gdel) # delete all entries that 
        r_starts=np.delete(r_starts,[i+1 for i in gdel])
        r_lens=r_ends-r_starts+1         # length of events recalculated


    # 
    print ("QUALITY CHECK ONM REVERSALS")
    for i,ilen in enumerate(r_lens):
        #print("event ",i)
        # smooth slope has to exceed 1 sigma in mean and event last more than a day to count:
        r_mean=np.mean(slopes_smth[r_starts[i]:r_ends[i]+1])
        r_before=np.mean(slopes_smth[r_starts[i]-24:r_starts[i]]) # day before
        r_after=np.mean(slopes_smth[r_ends[i]:r_ends[i]+24])      # day after
        if ilen>=pars["event_len"]:
            #print(slopes_smth[r_starts[i]:r_ends[i]+1])
            #print("length okay",ilen, "rmean ",r_mean, slopes_crit_lo)
            if r_mean<slopes_crit_lo:
                #print ("mean achieved ",r_mean)
                if r_before>slopes_crit_hi and r_after>slopes_crit_hi:
                    #print ("before and after positive")
                    if r_starts[i]>nextra and r_ends[i]<ntime-nextra:  # and any(slopes_mask_hi[r_starts[i]-48:r_starts[i]]) and any(slopes_mask_hi[r_ends[i]:r_ends[i]+48]):
                     # we have a good event, check the next is not straight after...
                        #print ("not near edge, event accepted!")
                        reverses.append({"start":r_starts[i],"end":r_ends[i],"len":ilen})

    nreverses=len(reverses)

    # here we construct the composite:

    fdate=date1.strftime("%Y%m%d")+"-"+date2.strftime("%Y%m%d")
    tagpdf="_"+region+"_anim"+fdate

    if not os.path.isdir(sdir["wrk"] + "utc9_time/sst/" + fdate):
        os.makedirs(sdir["wrk"] + "utc9_time/sst/" + fdate)
        
    if not os.path.isdir(sdir["wrk"] + "utc9_time/sst/slopes"):
        os.makedirs(sdir["wrk"] + "utc9_time/sst/slopes")

    # smooth the radiation fields to remove hte diurnal cycle
    nwin=72 # 24 
    sw_all=pd.Series(sw_all).rolling(window=nwin,min_periods=1,center=True).mean()
    sw_clr=pd.Series(sw_clr).rolling(window=nwin,min_periods=1,center=True).mean()

    #------------------------------------
    # dump to netcdf
    #------------------------------------
    var_slope = fileout.createVariable('slope', np.float64, 'time',)
    var_slope_np = fileout.createVariable('slope_np', np.float64, 'time',)
    var_slope.setncattr('units','m**-2')
    var_slope_np.setncattr('units','m**-2')
    var_slope[:]=slopes
    var_slope_np[:]=slopes_np
    fileout.close()


    # copy data to pandas dataframe - this is pretty YUKKY, it should be built from scratch.
    df["time"]=times
    df["slope"]=slopes
    df["year"]=years
    df["month"]=months
    df["day"]=days
    df["hour"]=hours
    df["mjo_amp"]=mjostrength
    df["mjo_phase"]=mjophase
    df["meansst"]=meansst
    df["minW"]=minW
    df["maxW"]=maxW
    df["meanW"]=meanW
    df["mld"]=mld
    df["minU"]=minu
    df["maxU"]=maxu
    df["meanU"]=avgu
    df["olr_clr"]=olr_clr
    df["radn_all"]=[sw_all[x]-olr_all[x] for x in range(len(olr_all))]
    df["precip"]=meanpp

    path=sdir["wrk"] + "utc9_time/sst/slopes/slope_"+fdate+".csv"
    df.to_csv(path,index=False)


    # Construct the composite
    complist=[
    {"name":"slope","units":"K %$^-1$","idx":0,"title":"slope"}
    ,{"name":"meansst","units":"K","idx":1,"title":"Sea Surface Temperature"}
    ,{"name":"minW","units":"kg m$^{-2}$","idx":2,"title":"Minimum Precipitable Water"}
#   ,{"name":"meanW","units":"kg m$^{-2}$","idx":2,"title":None}
#   ,{"name":"maxW","units":"kg m$^{-2}$","idx":2,"title":None}
    ,{"name":"meanU","units":"m s$^{-1}$","idx":3,"title":"Mean Zonal Wind"}
#    ,{"name":"maxU","units":"m s$^{-1}$","idx":3,"title":None}
    ,{"name":"olr_clr","units":"W m$^{-2}$","idx":4,"title":"Clear sky OLR"}
    ,{"name":"radn_all","units":"W m$^{-2}$","idx":5,"title":"Net TOA Radn"}
#    ,{"name":"precip","units":"mm hr$^{-1}$","idx":5,"title":"Precipitation"}
    ]
    npanelcomp=max([i["idx"] for i in complist])+1
    dfrev={}
    ncomp=2*nextra+ninterpol
    print ("ncomp",ncomp)
    print("WE FOUND THIS MANY EVENTS: ",nreverses)

    for var in complist:
        dfrev[var["name"]]=np.full([ncomp,nreverses],np.nan)
        for irev,rev in enumerate(reverses):
            xp=np.array(range(rev["len"]))
            xi=(rev["len"]-1)*np.linspace(0,1,ninterpol)
            print ( df[var["name"]][rev["start"]:rev["end"]+1])
            newdata=np.interp(xi,xp,df[var["name"]][rev["start"]:rev["end"]+1])

            # paste on start and end, N=48 hourly datapoints.
            dfrev[var["name"]][:,irev]=np.concatenate([
                       df[var["name"]].to_numpy()[rev["start"]-nextra:rev["start"]],
                       newdata,
                       df[var["name"]].to_numpy()[rev["end"]+1:rev["end"]+nextra+1]])

    # average the events 
    for var in complist:
        dfrev[var["name"]]=np.nanmean(dfrev[var["name"]],axis=1)
        #print(dfrev[var["name"]].shape)


    #
    # Composite plot
    #
    fig,(ax)=plt.subplots(npanelcomp,1,sharex=True,figsize=(10,6))
    for ivar,var in enumerate(complist):
        idx=var["idx"]
        ax[idx].plot(range(ncomp),dfrev[var["name"]])
        if var["title"]:
            ax[idx].set_title("("+ascii_lowercase[idx]+") "+var["title"],fontsize=10)
        if ivar==len(complist)-1:
            ax[idx].set_xlabel("Days")
        ax[idx].set_ylabel(var["units"])
        ax[idx].axvline(x=nextra,linewidth=1.5, linestyle='--', color='red')
        ax[idx].axvline(x=ncomp-nextra,linewidth=1.5, linestyle='--', color='red')
        xticks=list(range(0,nextra+1,24))+[ncomp/2]+list(range(ncomp-nextra,ncomp+1,24))
        xlabs=[str(i) for i in range(-pars["nextrad"],0)]+["start","Reversal","end"]+[str(i) for i in range(1,pars["nextrad"]+1)]
        ax[idx].set_xticks(xticks) # choose which x locations to have ticks
        ax[idx].set_xticklabels(xlabs) # set the labels to display at those ticks
        
        #ax[idx].text(nextra/2,np.min(dfrev[var["name"]]),"days before event",horizontalalignment="center")
        #ax[idx].text(ncomp-nextra/2,np.min(dfrev[var["name"]]),"days after event",horizontalalignment="center")
    fig.tight_layout()
    plt.savefig(sdir["plt"]+'composite_gpm'+tagpdf+plt_type, dpi=300)
    plt.close()
    
    print ("here ",sdir["plt"]+'composite_gpm'+tagpdf+plt_type)
    #
    # this is the hovmuller sst as function of PW percentile
    #

    lwd=0.8
    panwd=1.0
    ip=0 # plot index

    for var in ['sst']:
        # HARDWIRED YUK!
        fig,(ax)=plt.subplots(1,8,figsize=(15,10),gridspec_kw={'width_ratios': [5]+[panwd]*7},sharey=True)
        fig.subplots_adjust(left=0.15)
        #plt.rc('font', family='serif')
        #ax.margins(0.2)
        img=ax[ip].pcolormesh(hovx,times,hov[var][1:ntime+1],cmap='seismic',rasterized=True,vmin=sstamin,vmax=sstamax)
        preciplevs=[0.5,5]
        ax[ip].contour(hovx,times,hov['precip'][1:,:],levels=preciplevs,colors=['black'],linewidths=[1.0,2.0])
        ax[ip].contourf(hovx,times,hov['precip'][1:,:],levels=preciplevs,colors='none',hatches=['....','////'])
        ax[ip].set_xlabel('TPW percentile', fontsize=12)
        ax[ip].set_ylabel('Date (UTC+09:00)',fontsize=12)

        ax[ip].set_yticks(ytick)
        ax[ip].set_yticklabels(yticklabels,fontsize=10)
        cbar=plt.colorbar(img,ax=ax[0])
        cbar.ax.tick_params(labelsize=9) 
        #cbar.set_label('SST Anomaly ($^o$C)',fontsize=12,rotation=-90,y=1)
        ax[ip].set_title('SST Anomaly ($^o$C)',loc='right')

        ip+=1

        # Add time series?
        ax[ip].plot(sst_hires_min,times, color='blue', linewidth=lwd)
        ax[ip].plot(sst_hires_mean,times, color='lightgreen', linewidth=lwd)
        ax[ip].plot(sst_hires_max,times, color='red', linewidth=lwd)
        ax[ip].set_xlabel('SST H8 ($^o$C)', fontsize=8)
        add_rev(ax[ip],np.min(sst_hires_min),np.max(sst_hires_max),reverses)
        #ax[ip].set_xlim(29,31)
        ax[ip].tick_params(axis='x', labelsize=5)
        ax[ip].legend(['min', 'mean', 'max'], bbox_to_anchor=(0, 1), loc='lower left', fontsize=4)
        
        ip+=1

        ax[ip].plot(slopes,times,linewidth=lwd)
        ax[ip].plot(slopes_smth,times,color="red",linewidth=0.5)
        ax[ip].plot(slopes_np_smth,times,color="purple",linewidth=0.5)

        add_rev(ax[ip],np.min(slopes),np.max(slopes),reverses)
        ax[ip].axvline(x=0, linewidth=lwd, linestyle='--', color='0.8')
        ax[ip].axvline(x=slopes_av, linewidth=lwd, color='orange')
        ax[ip].axvline(x=slopes_crit_lo, linewidth=lwd, linestyle=':', color='orange')
        ax[ip].set_xlabel('slope', fontsize=8)
        ax[ip].tick_params(axis='x', labelsize=3)
        
        ip+=1

        ax[ip].plot(meanpp,times)
        ax[ip].set_xlabel('mean Precip (mm/hr)', fontsize=7)
        add_rev(ax[ip],np.min(meanpp),np.max(meanpp),reverses)
        ax[ip].tick_params(axis='x', labelsize=5)
  
        ip+=1

        # AMT: I replaced sigma TPW with MLD
        #ax[ip].plot(mld,times)
        #ax[ip].set_xlabel('MLD (m)', fontsize=8)
        #add_rev(ax[ip],np.min(mld),np.max(mld),reverses)
        #ax[ip].tick_params(axis='x', labelsize=5)
        #ip+=1

        ax[ip].plot(df["minW"],times, color='blue', linewidth=lwd)
        ax[ip].plot(df["meanW"],times, color='lightgreen', linewidth=lwd)
        ax[ip].plot(df["maxW"],times, color='red', linewidth=lwd)
        ax[ip].axvline(x=np.nanmean(minW), linewidth=lwd, linestyle='--', color='0.8')
        ax[ip].axvline(x=np.nanmean(maxW), linewidth=lwd, linestyle='--', color='0.8')
        add_rev(ax[ip],np.min(df["minW"]),np.max(df["maxW"]),reverses)
        ax[ip].set_xlabel('TPW (mm)', fontsize=8)
        ax[ip].tick_params(axis='x', labelsize=5)
        ax[ip].legend(['min', 'mean', 'max'], bbox_to_anchor=(0, 1), loc='lower left', fontsize=4)
        ip+=1

# =============================================================================
#         zip_object = zip(minW, maxW)
#         for list1_i, list2_i in zip_object:
#             diff.append(list1_i-list2_i)
#         res5=ax[4].plot(np.abs(diff),times)
#         ax[4].set_xlabel('maxW-minW', fontsize=8)
#         ax[4].tick_params(axis='x', labelsize=5)
# =============================================================================
        
        ax[ip].plot(minu,times, color='blue', linewidth=lwd)
        ax[ip].plot(avgu,times, color='lightgreen', linewidth=lwd)
        ax[ip].plot(maxu,times, color='red', linewidth=lwd)
        ax[ip].axvline(x=0, linewidth=lwd, linestyle='--', color='0.8')
        add_rev(ax[ip],np.min(minu),np.max(maxu),reverses)
        ax[ip].set_xlabel('10m U(m s$^{-1}$)', fontsize=8)
        ax[ip].tick_params(axis='x', labelsize=5)
        ax[ip].legend(['min', 'mean', 'max'], bbox_to_anchor=(0, 1), loc='lower left', fontsize=4)
        ip+=1

        #im1=ax[7].scatter(mjostrength,times, c=mjophase, cmap='Accent', s=1, vmin=0.5, vmax=8.5)
        #ax[7].axvline(x=1, linewidth=lwd, linestyle='--', color='0.8')
        #ax[7].set_xlabel('OMI', fontsize=8)
        #ax[7].tick_params(axis='x', labelsize=4)
        #cbar2=plt.colorbar(im1, ax=ax[7])
        #cbar2.ax.tick_params(labelsize=5)

 
        ax[ip].plot(olr_all,times,color="red",linewidth=lwd)
        ax[ip].plot(olr_clr,times,color="blue",linewidth=lwd)
        #ax[ip].plot(sw_all,times,color="blue",linewidth=lwd)
        #ax[ip].plot(sw_clr,times,color="blue",linewidth=lwd)
        #ax[ip].legend(["lw all","lw clr","sw all","sw clr"],bbox_to_anchor=(0, 1), loc='lower left', fontsize=4)
#        ax[ip].legend(["lw all","sw all"],bbox_to_anchor=(0, 1), loc='lower left', fontsize=4)
        add_rev(ax[ip],np.min(olr_clr),np.max(olr_clr),reverses)
        ax[ip].set_xlabel('Flux (W m$^{-2}$)', fontsize=8)
        ip+=1
  
        # total net radn 
        tot_all=[sw_all[x]-olr_all[x] for x in range(len(olr_all))]
        tot_clr=[sw_clr[x]-olr_clr[x] for x in range(len(olr_clr))]

        ax[ip].plot(tot_all,times,color="red",linewidth=lwd)
        ax[ip].set_xlim([40,160])
#        ax[ip].plot(tot_clr,times,color="red",linewidth=lwd)
        ax[ip].legend(["all"],bbox_to_anchor=(0, 1), loc='lower left', fontsize=4)
        ax[ip].set_xlabel('LW+SW (W m$^{-2}$)', fontsize=8)
        ip+=1

        #plt.savefig(mdir.dods+'sst/hov'+var+'_gpm'+tagpdf+plt_type)
        fig.tight_layout()
        plt.savefig(sdir["plt"]+'hov'+var+'_gpm'+tagpdf+plt_type, dpi=300)
        plt.close()

    # line plots of precip as function of vars
    for var in vars:
        pav=np.divide(ptot[var], ntot[var],out=np.zeros_like(ptot[var])-99, where=ntot[var]!=0)
        pav=np.nan_to_num(pav)
        plt.plot(bins[var][0:nbins[var]],pav)
        plt.ylim(0,np.max(pav))
        plt.ylabel('Precip (mm/hr)')
        plt.title('Western Pacific')
        plt.xlabel(vartitle[var])
        #plt.savefig(mdir.dods+'sst/'+var+'_all_gpm'+tagpdf+plt_type)
        plt.savefig(sdir["plt"]+var+'_all_gpm'+tagpdf+plt_type)
        plt.close()

        pav=np.divide(ppostot[var],npostot[var],out=np.zeros_like(ptot[var])-99, where=ntot[var]!=0)
        pav=np.nan_to_num(pav)
        plt.plot(bins[var][0:nbins[var]],pav)
        plt.ylim(0,np.max(pav))
        plt.ylabel('Precip when raining (mm hr$^{-1}$)')
        plt.title('Western Pacific')
        plt.xlabel(vartitle[var])
        #plt.savefig(mdir.dods+'sst/'+var+'_pos_gpm'+tagpdf+plt_type)
        plt.savefig(sdir["plt"]+var+'_pos_gpm'+tagpdf+plt_type)
        plt.close()

        pav=np.divide(npostot[var],ntot[var],out=np.zeros_like(ptot[var])-99, where=ntot[var]!=0)
        pav=np.nan_to_num(pav)
        plt.plot(bins[var][0:nbins[var]],pav)
        plt.ylim(0,np.max(pav))
        plt.ylabel('Precip Frequency (mm hr$^{-1}$)')
        plt.title('Western Pacific')
        plt.xlabel(vartitle[var])
        #plt.savefig(mdir.dods+'sst/'+var+'_freq_gpm'+tagpdf+plt_type)
        plt.savefig(sdir["plt"]+var+'_freq_gpm'+tagpdf+plt_type)
        plt.close()

    #
    # block plots of PW real vals, this is slow so last
    #
    for var in ['sst']:
        fig,(ax)=plt.subplots(1,2,gridspec_kw={'width_ratios': [4, 1]},sharey=True,squeeze=False)
        #fig,(ax)=plt.subplots(1,2,squeeze=False)
        ax[0,0].set_ylim(0,ntime)
        fig.subplots_adjust(left=0.2)
        ax[0,0].set_xlim(20,80)

        cmap=plt.cm.seismic  
        cols=cmap((hov['sst'][0:ntime,:]-sstamin)/(sstamax-sstamin))

        start_time = timeit.default_timer()
        
        for itime in range(1,ntime):
#           for ihov in range(nhov):
#            lpatch=[patches.Rectangle((hovmin['tpw'][itime,ihov],itime),hovmax['tpw'][itime,ihov]-hovmin['tpw'][itime,ihov],1.0,facecolor=cols[itime,ihov]) for itime in range(1,ntime)]
            lpatch=[patches.Rectangle((hovmin['tpw'][itime,ihov],itime),hovmax['tpw'][itime,ihov]-hovmin['tpw'][itime,ihov],1.0,facecolor=cols[itime,ihov]) for ihov in range(nhov)]
            collection = collections.PatchCollection(lpatch,match_original=True)
            ax[0,0].add_collection(collection) 

        ax[0,0].set_xlabel('PW (kg m$^{-2}$)', fontsize=13)
        ax[0,0].set_ylabel('Date ',fontsize=15)
        ax[0,0].set_yticks(ytick)
        ax[0,0].set_yticklabels(yticklabels,fontsize=8)
        cbar=plt.colorbar(img,ax=ax[0,0])
        #cbar.set_label('SST Anomaly ($^o$C)',fontsize=15)
        ax[0,0].set_title('SST Anomaly ($^o$C)',loc='right')

        # Add time series?
        res=ax[0,1].plot(meansst,times)
        ax[0,1].set_xlabel('SST ($^o$C)', fontsize=12)
        ax[0,1].set_xlim(sstmin,sstmax)
        
        #plt.savefig(mdir.dods+"sst/blockplot_distall_"+var+'_gpm'+tagpdf+plt_type)
        plt.savefig(sdir["plt"]+"blockplot_distall_"+var+'_gpm'+tagpdf+plt_type)
        plt.close()
        elapsed = timeit.default_timer() - start_time
        print("******** block time ",elapsed)


def main(argv):
    # Get command line options
    global pars, sdir
    pars = setup.defaults() 
    pars = getargs(pars)

    sdir=locs(pars)
    print (sdir["wrk"])
    print('starting')
    startdates = pd.date_range(pars["date1"],pars["date2"])
    for date in startdates:
       print (date)
       oneplot(date,pars)

if __name__ == "__main__":
    main(sys.argv[1:])





