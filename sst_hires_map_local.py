import os
os.environ["PROJ_LIB"] = r'C:\Users\Michie\anaconda3\Library\share\basemap'
from netCDF4 import Dataset, num2date, date2index
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
import dateutil.parser
from datetime import datetime
from matplotlib import colors
import gc
from scipy.spatial import KDTree


llon = 135
rlon = 145

# =============================================================================
# llat = 5
# ulat = 15
# =============================================================================

llat = 2
ulat = 9

#data_dir = "/media/mde_vera/Seagate Backup Plus Drive/w_"+str(llat)+"_"+str(ulat)+"_data"
#data_dir = "D:\w_"+str(llat)+"_"+str(ulat)+"_data"
data_dir = "D:\\"+str(llat)+"_"+str(ulat)+"_data"

folder_name = "filled_sst_hires_uv_maps"
var = "sst_hires"
var2 = "cld"
var3 = "mim"

utc=9

dates=[("20170201", "20170228"), ("20170301", "20170331"), ("20170401", "20170430")]


#("20170101", "20170131"), ("20170201", "20170228"), ("20170301", "20170331"), ("20170401", "20170430")
#("20170501", "20170531") 
#("20170601", "20170630"),("20170701", "20170731"), ("20170801", "20170831"), 
#("20170901", "20170930"),("20171001", "20171031"), ("20171101", "20171130")
#("20171201", "20171231")


for (date1, date2) in dates:
    daterange = pd.date_range(date1, date2)
    fdate=date1+"-"+date2
    #print(daterange)
    
    if not os.path.isdir(data_dir + "/utc9_time/" + folder_name):
        os.makedirs(data_dir + "/utc9_time/" + folder_name)
    
    if not os.path.isdir(data_dir + "/utc9_time/" + folder_name + "/" + fdate):
        os.makedirs(data_dir + "/utc9_time/" + folder_name+ "/" + fdate)
        
    for idate,sdate in enumerate(daterange):
        date=sdate.strftime("%Y%m%d")
        yy=sdate.strftime("%Y")
        mm=sdate.strftime("%m")
        dd=sdate.strftime("%d")
        
        locdatestart=sdate+pd.DateOffset(hours=-utc)
        locdateend=locdatestart+pd.DateOffset(hours=23)
        locdaterange = pd.date_range(locdatestart,locdateend, freq='H')
        print(locdaterange)
    
        fig, axs = plt.subplots(4, 6, figsize=(60,40))
        
        for ndate,ldate in enumerate(locdaterange):
            fulldate=ldate.strftime("%Y%m%d%H")
            lyy=ldate.strftime("%Y")
            lmm=ldate.strftime("%m")
            ldd=ldate.strftime("%d")
            lhh=ldate.strftime("%H")
            print(fulldate)
            
            dfilename = data_dir+"/"+var+"_wpac_"+str(llat)+"-"+str(ulat)+"N_"+fulldate+".nc"
            try:
                dset = Dataset(dfilename)
                dvar = dset.variables["sea_surface_temperature"][:]
                #dvar = dset.variables["l2p_flags"][:]
            except:
                print("Error opening file ", dfilename)
    
            lon = dset.variables["lon"][:]
            lat = dset.variables["lat"][:]
    
            dfilename2 = data_dir+"/"+var2+"_wpac_"+str(llat)+"-"+str(ulat)+"N_"+fulldate+".nc"
            try:
                dset2 = Dataset(dfilename2)
            except:
                print("Error opening file ", dfilename2)
    
            dvar2 = dset2.variables["tope"][:]
    
     
            masked_dvar = np.ma.array(dvar,mask=np.isnan(dvar))
            masked_dvar2 = np.ma.array(dvar2,mask=np.isnan(dvar2))
            
            missing_index=np.logical_and(masked_dvar.mask, masked_dvar2.mask)
            #print(missing_index[0])

            x,y=np.mgrid[0:dvar[0].shape[0],0:dvar[0].shape[1]]
            #print(x)
            #print(y)
            
            xygood = np.array((x[~missing_index[0]],y[~missing_index[0]])).T
            #print(xygood)
            xybad = np.array((x[missing_index[0]],y[missing_index[0]])).T
            #print(xybad)
            
            #print(KDTree(xygood).query(xybad)[1])
            
            dvar[missing_index]=dvar[~missing_index][KDTree(xygood).query(xybad)[1]]
            #print(dvar)

    
            
    # =============================================================================
    #         dfilename3 = data_dir+"/"+var3+"_wpac_"+str(llat)+"-"+str(ulat)+"N_"+date+str(hr).zfill(2)+".nc"
    #         try:
    #             dset3 = Dataset(dfilename3)
    #             dvar3 = dset3.variables["tpwGrid"][:]
    #         except:
    #             print("Error opening file ", dfilename3)
    # =============================================================================
    
            
            #print(dset.variables["l2p_flags"][:].getncattr("flag_masks"))
            #print(dvar)
            #binary_repr_v = np.vectorize(np.binary_repr)
            #print(binary_repr_v(dvar, 15))
            #print(dvar)
            
            #u wind data
            uwindfile=data_dir+'/10m_u'+"_wpac_"+str(llat)+"-"+str(ulat)+"N_"+lyy+".nc"
            try:
                duwind=Dataset(uwindfile) 
            except:
                print("file error:",uwindfile)
            
            timenow=ldate.strftime("%Y-%m-%dT%H:00:00")
            dt = dateutil.parser.parse(timenow)
            
            u_all_times=duwind.variables["time"]
            u_dt_idx = date2index(dt, u_all_times)
            
            uval=duwind.variables['u10'][u_dt_idx]
            #print(uval.shape)
            #print(uval[0:100:10,0:100:10])
            
            #v wind data
            vwindfile=data_dir+'/10m_v'+"_wpac_"+str(llat)+"-"+str(ulat)+"N_"+lyy+".nc"
            try:
                dvwind=Dataset(vwindfile) 
            except:
                print("file error:",vwindfile)
            
            v_all_times=dvwind.variables["time"]
            v_dt_idx = date2index(dt, v_all_times)
            
            vval=dvwind.variables['v10'][v_dt_idx]
    
    
            mdate=ldate+pd.DateOffset(hours=utc)
            #mapdate=mdate.strftime("%Y%m%d%H")
            hh=mdate.strftime("%H")
            
            row = int(int(hh)/6)
            col = int(hh)%6
            axs[row][col].set_title(str(hh).zfill(2), fontsize = 45)
            m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=rlon, resolution='i', ax=axs[row][col])
            m.drawcoastlines()
            m.drawparallels(np.arange(llat,ulat,2.), labels=[1,0,0,0], fontsize = 25)
            m.drawmeridians(np.arange(llon,rlon+1,2.), labels=[0,0,0,1], fontsize = 25)
    
            lon, lat = np.meshgrid(lon, lat)
            x, y = m(lon, lat)
            #cbticks = np.arange(300, 305.5, 0.5) 
    
           
            #temp = m.scatter(x, y, c=dvar[0,:,:], cmap = "RdYlBu_r")
            temp = m.contourf(x, y, dvar[0,:,:], levels=np.arange(301, 305.5, 0.5), alpha = 0.9, cmap = "RdYlBu_r", extend="both")
            
            height = m.contourf(x, y, dvar2[0,:,:]/100, levels=np.arange(0, 21, 10), colors=['lightgrey', 'grey'], alpha = 0.3) 
            #height = m.contour(x, y, dvar2[0,:,:]/100, vmin=10, colors='m', alpha = 0.7)
            

            #try:
            #    water = m.contourf(x, y, dvar3, levels=np.arange(20,90,5), hatches=['---', '-', '//', '/', './', '\\ \\', '\\', '\.', '*', '-|', '|', '.|', '.-', '...'], colors='none', extend="both")
            #except:
            #    print("Error water plot ", dfilename3)
            #divnorm = colors.TwoSlopeNorm(vmin=-10., vcenter=0, vmax=10)
            #uvwind = m.quiver(x, y, uval, vval, uval, cmap = "seismic", norm=divnorm)
            uvwind = m.quiver(x[0:70:4,0:100:9], y[0:70:4,0:100:9], uval[0:70:4,0:100:9], vval[0:70:4,0:100:9])
            
            #ax.clabel(height, inline=True, fontsize=10)
            #cs = map.contour(xx, yy, data, range(400, 1500, 100), cmap = plt.cm.cubehelix)
            #plt.clabel(height, inline=True, fmt='%1.0f', fontsize=12, colors='k')
            #cb = m.colorbar(temp, "right", size="5%", pad="2%")
            #plt.clim(270,300)
            #cb.set_label('Temperature (K)')
    
        fig.suptitle("Himawari 8 SST: " + date, fontsize = 70)
        cax = fig.add_axes([0.925,0.1,0.015,0.8])
        cb = fig.colorbar(temp, cax = cax)
        cb.ax.tick_params(labelsize = 40)
        cb.ax.get_yaxis().labelpad = 60
        cb.ax.set_ylabel('Temperature (K)', rotation=270, fontsize = 50)
        #plt.show()
        
    # =============================================================================
    #     cax2 = fig.add_axes([0.075,0.1,0.015,0.8])
    #     cb2 = fig.colorbar(water, cax = cax2)
    #     cb2.ax.tick_params(labelsize = 40)
    #     cb2.ax.get_yaxis().labelpad = 60
    #     cb2.ax.yaxis.set_ticks_position('left')
    #     cb2.ax.yaxis.set_label_position('left')
    #     cb2.ax.set_ylabel('Total Precipitable Water (mm)', rotation=90, fontsize = 50)
    # =============================================================================
    
        plt.savefig(data_dir+"/utc9_time/"+folder_name+"/"+fdate+"/filled_"+var+"_cth_"+str(llat)+"-"+str(ulat)+"N_"+date+".png")
        plt.cla()
        plt.clf()
        plt.close("all")
        plt.close(fig)
        gc.collect()

