import xarray as xr
import matplotlib.pyplot as plt
import pdb

#Paths!
reals = '/home/tompkins-archive/acasallas/Real1_run/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Reals/'
scra = '/home/netapp-clima/scratch/acasallas/wrf/Real_run/WRF/'

####################################################################
### Read data
bot = '0'
time = 55

cont = 'control'
print('Reading control')
ds_c = xr.open_dataset(reals+cont+'/test_'+cont+'.nc')
ps = ds_c.PSFC[time,:,:]
znu =  ds_c.ZNU[0,:]*ps[0,0]/100.

test = '20_mar'
print('Reading 20 mar')
ds_org = xr.open_dataset(reals+cont+'/test_'+cont+'.nc') 
print('Reading 12 apr')
ds_12 = xr.open_dataset(scra+'/test_12_apr_file02.nc')

print('Selecting variables')
con_u = ds_c.U
org_u = ds_org.U
rev_u = ds_12.U
org_v = ds_org.V
rev_v = ds_12.V

del(ds_c);del(ds_org);del(ds_12)

### Select days
print('Selecting reversal days')
# Reversals
sta = [2*24, 10*24, 27.5*24, 36*24, 41*24, 47*24]
end = [4.7*24, 11.5*24, 30.2*24, 37*24, 44*24, 53*24]

#r1 = con_u[int(sta[0]):int(end[0]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#r2 = con_u[int(sta[1]):int(end[1]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#r3 = con_u[int(sta[2]):int(end[2]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#r4 = con_u[int(sta[3]):int(end[3]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#r5 = con_u[int(sta[4]):int(end[4]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#r6 = con_u[int(sta[5]):int(end[5]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])

#rt = (r1+r2+r3+r4+r5+r6)/6
#del(r1);del(r2);del(r3);del(r4);del(r5);del(r6)

uo = org_u[int(20*24):int(25*24),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
ur = rev_u[int(10*24):int(15*24),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
vo = org_v[int(20*24):int(25*24),:,:,:].mean(dim=['Time','south_north_stag','west_east'])
vr = rev_v[int(10*24):int(15*24),:,:,:].mean(dim=['Time','south_north_stag','west_east'])
### Organize days
print('Selecting organize days')
sta = [6*24,12*24,16*24,34*24,45*24]
end = [7*24,14*24,25*24,35*24,46*24]

#o1 = con_u[int(sta[0]):int(end[0]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#o2 = con_u[int(sta[1]):int(end[1]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#o3 = con_u[int(sta[2]):int(end[2]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#o4 = con_u[int(sta[3]):int(end[3]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])
#o5 = con_u[int(sta[4]):int(end[4]),:,:,:].mean(dim=['Time','south_north','west_east_stag'])

#ot = (o1+o2+o3+o4+o5)/5
#del(o1);del(o2);del(o3);del(o4);del(o5)

plt.figure(figsize=(6,9))
#plt.plot(rt,znu, label = 'Reversals', color = 'b')
#plt.plot(ot,znu, label = 'Organize', color = 'purple')
plt.plot(uo,znu, label = 'Organize run', color = 'purple')
plt.plot(ur,znu, label = 'Reversal run', color = 'blue')
plt.ylim(1020,100)
plt.xlim(-12,0)
plt.legend()
plt.ylabel('Pressure (hPa)')
plt.xlabel('Zonal wind (m s$^{-1}$)')
plt.savefig(med+'Zonal_shear_workshop.jpg', bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(6,9))
#plt.plot(rt,znu, label = 'Reversals', color = 'b')
#plt.plot(ot,znu, label = 'Organize', color = 'purple')
plt.plot(vo,znu, label = 'Organize run', color = 'purple')
plt.plot(vr,znu, label = 'Reversal run', color = 'blue')
plt.ylim(1020,100)
plt.legend()
plt.ylabel('Pressure (hPa)')
plt.xlabel('Meridional wind (m s$^{-1}$)')
plt.savefig(med+'Meridional_shear_workshop.jpg', bbox_inches='tight')
plt.show()
plt.close()

pdb.set_trace()

