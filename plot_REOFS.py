import pandas as pd
import xarray as xr
import numpy as np
import pdb
from shapely.geometry.polygon import LinearRing
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as GridSpec
from pylab import *
import os
from os.path import exists

scra = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/'
user = '/home/netapp-clima/users/acasallas/ERA5_real/'

ds_sst = xr.open_dataset(user+'REOFs_SST_all.nc')
ds_tcwv = xr.open_dataset(user+'REOFs_TCWV_all.nc')

fig = plt.figure(figsize=(14,8))
gs = GridSpec(5,3,left = 0.1, right = 0.93, hspace=0.25, wspace=0.1, top = 0.92, bottom = 0.08,height_ratios=[1,0.1,0.1,1,0.1])

#######################
levtcw = [0,0.005,0.01,0.015,0.02,0.025,0.03,0.035]
ax = plt.subplot(gs[0],projection= ccrs.PlateCarree(central_longitude=180))
im = plt.contourf(ds_tcwv['lon'],ds_tcwv['lat'],ds_tcwv.EOFs[:,:,0], cmap = 'winter_r', transform = ccrs.PlateCarree(), levels=levtcw, extend='both')
gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1, color='k', alpha = 0.5, linestyle = ':')
gl.xlocator = mticker.FixedLocator([137.5, 140, 142.5]);gl.xlabels_bottom = False
gl.ylocator = mticker.FixedLocator([2, 3.75, 5.5, 7.25, 9]);gl.ylabels_right = False
plt.title('(a) TCWV REOF mode 1', fontweight = 'bold', loc = 'left')

#######################
ax = plt.subplot(gs[1],projection= ccrs.PlateCarree(central_longitude=180))
im = plt.contourf(ds_tcwv['lon'],ds_tcwv['lat'],ds_tcwv.EOFs[:,:,1], cmap = 'winter_r', transform = ccrs.PlateCarree(), levels=levtcw, extend='both')
gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1, color='k', alpha = 0.5, linestyle = ':')
gl.xlocator = mticker.FixedLocator([137.5, 140, 142.5]);gl.xlabels_bottom = False
gl.ylocator = mticker.FixedLocator([2, 3.75, 5.5, 7.25, 9]); gl.ylabels_right = False
gl.ylabels_left = False
plt.title('(b) TCWV REOF mode 2', fontweight = 'bold', loc = 'left')

#######################
levtcw = [-0.03,-0.025,-0.02,-0.015,-0.01,-0.005,0]
ax = plt.subplot(gs[2],projection= ccrs.PlateCarree(central_longitude=180))
im = plt.contourf(ds_tcwv['lon'],ds_tcwv['lat'],ds_tcwv.EOFs[:,:,2], cmap = 'winter_r', transform = ccrs.PlateCarree(), levels=levtcw, extend='both')
gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1, color='k', alpha = 0.5, linestyle = ':')
gl.xlocator = mticker.FixedLocator([137.5, 140, 142.5]); gl.xlabels_bottom = False
gl.ylocator = mticker.FixedLocator([2, 3.75, 5.5, 7.25, 9]); gl.ylabels_left = False
plt.title('(c) TCWV REOF mode 3', fontweight = 'bold', loc = 'left')

ax =  plt.subplot(gs[1,:])
plt.colorbar(im, cax=ax, orientation = 'horizontal')

#######################
sstlev = [-5e-5,0.0,0.00005,0.00010,0.00015,0.0002]
ax = plt.subplot(gs[9],projection= ccrs.PlateCarree(central_longitude=180))
im1 = plt.contourf(ds_sst['lon'],ds_sst['lat'],ds_sst.EOFs[:,:,0], cmap = 'autumn_r', transform = ccrs.PlateCarree(), levels=sstlev, extend='both')
gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1, color='k', alpha = 0.5, linestyle = ':')
gl.xlocator = mticker.FixedLocator([137.5, 140, 142.5]); gl.xlabels_top = False
gl.ylocator = mticker.FixedLocator([2, 3.75, 5.5, 7.25, 9]); gl.ylabels_right = False
plt.title('(d) SST REOF mode 1', fontweight = 'bold', loc = 'left')

#######################
ax = plt.subplot(gs[10],projection= ccrs.PlateCarree(central_longitude=180))
im1 = plt.contourf(ds_sst['lon'],ds_sst['lat'],ds_sst.EOFs[:,:,1], cmap = 'autumn_r', transform = ccrs.PlateCarree(), levels=sstlev, extend='both')
gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1, color='k', alpha = 0.5, linestyle = ':')
gl.xlocator = mticker.FixedLocator([137.5, 140, 142.5]); gl.xlabels_top = False 
gl.ylocator = mticker.FixedLocator([2, 3.75, 5.5, 7.25, 9]); gl.ylabels_right = False
gl.ylabels_left = False
plt.title('(e) SST REOF mode 2', fontweight = 'bold', loc = 'left')

ax = plt.subplot(gs[4,0:2])
plt.colorbar(im1, cax=ax, orientation = 'horizontal')

#######################
sstlev = [-0.0001,-0.00005,0.0,0.00005,0.0001]
ax = plt.subplot(gs[11],projection= ccrs.PlateCarree(central_longitude=180))
im1 = plt.contourf(ds_sst['lon'],ds_sst['lat'],ds_sst.EOFs[:,:,3], cmap = 'autumn_r', transform = ccrs.PlateCarree(), levels=sstlev, extend='both')
gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, linewidth = 1, color='k', alpha = 0.5, linestyle = ':')
gl.xlocator = mticker.FixedLocator([137.5, 140, 142.5]); gl.xlabels_top = False
gl.ylocator = mticker.FixedLocator([2, 3.75, 5.5, 7.25, 9]); gl.ylabels_left = False
plt.title('(f) SST REOF mode 3', fontweight = 'bold', loc = 'left')

ax = plt.subplot(gs[4,2])
plt.colorbar(im1, cax=ax, orientation = 'horizontal')
plt.savefig(scra+'REOFs_SST_TCWV.jpg', dpi = 300, bbox_inches = 'tight')
plt.savefig(scra+'REOFs_SST_TCWV.pdf', dpi = 300, bbox_inches = 'tight')
plt.show()
