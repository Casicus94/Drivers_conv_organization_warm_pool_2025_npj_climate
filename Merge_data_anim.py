import xarray as xr
import numpy as np
import pdb

def create_dict(varis):
    dict_list = {}
    for var in varis:
        dict_list[var] = [] 

    return(dict_list)

lat = '3-10'
lon = '147-157'
scra = '/home/netapp-clima/scratch/acasallas/Composites_'+lat+'/'
mdir = '/home/tompkins-archive2/acasallas/ERA5_real/Area_'+lat+'_'+lon+'/Anim_data/'

varis = ['SST','TCWV','OLR','mu','mv']
data_rev = create_dict(varis)
data_org = create_dict(varis)
varname = ['sst', 'tcwv', 'ttr', 'p71.162', 'p72.162']
tipo = '' #''


if tipo == '':
    seas = []
    #9.00001
    for i in np.arange(-9,9.00001,1/24):
        seas.append(str(i))
else:
    seas = ['9','8.5','8','7.5','7','6','5','4.5','4','3.5','3','2','1','0','-1','-2','-3']

for i,var in enumerate(varis):
    print('Merging '+var)
    try:
        ds = xr.open_dataset(mdir+var+'_full_rev_composite_'+lat+tipo+'.nc')
    except:
        for sea in seas:
            ds = xr.open_dataset(scra+var+'_rev_'+sea+'_'+lat+'.nc')
            data_rev[var].append(ds[varname[i]])
            ds = xr.open_dataset(scra+var+'_org_'+sea+'_'+lat+'.nc')
            data_org[var].append(ds[varname[i]])
        ds_rev = xr.concat(data_rev[var], dim='time')
        if var == 'TCWV' or var == 'OLR':
            ds_rev = ds_rev.drop_vars('season')
        ds_rev.to_netcdf(mdir+var+'_full_rev_composite_'+lat+tipo+'.nc')
        ds_org = xr.concat(data_org[var], dim='time')
        if var == 'TCWV' or var == 'OLR':
            ds_org = ds_org.drop_vars('season')
        ds_org.to_netcdf(mdir+var+'_full_org_composite_'+lat+tipo+'.nc')

