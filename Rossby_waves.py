import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from datetime import datetime, timedelta
from scipy.signal import savgol_filter
import pdb

print('Starting Wave Script!, Libraries loaded')
# Some functions from AMS
def get_minima_indicies(array):
    return (np.diff(np.sign(np.diff(array))) > 0).nonzero()[0] + 1

def get_maxima_indicies(array):
    return (np.diff(np.sign(np.diff(array))) < 0).nonzero()[0] + 1

def get_increments(previous_extrema, current_extrema):
    # Make sure these are cast as lists
    previous_extrema = list(previous_extrema)
    current_extrema = list(current_extrema)
    # Adjust for overflow
    current_extrema = cycle_adjust(previous_extrema, current_extrema)
    increments = []
    for lon in current_extrema:
        # Find closest extrema in previous, and assume continuation
        prev_lon = min(previous_extrema, key=lambda x:abs(x-lon))
        increments.append(lon - prev_lon)
    return np.array(increments)

def cycle_adjust(previous_extrema, current_extrema):
    # Bump up an overflowed element
    if min(current_extrema) < min(previous_extrema)/2:
        # If the the new lowest element is closer to zero than it is to the previous lowest,
        # assume it is an overflow
        overflowed = current_extrema.pop(0)
        current_extrema.append(overflowed + 360.)
    return current_extrema

def reject_outliers(data):
    d = data - np.median(data)
    iqr = np.percentile(d, 75) - np.percentile(d, 25)
    iqr = iqr if iqr > 0 else 1.
    return data[np.abs(d) < 1.5 * iqr]

def fraction_of_day(td):
    return td / np.timedelta64(1, 'D')

##Paths
scr = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/ERA5/'
user = '/home/netapp-clima/users/acasallas/ERA5_real/'
print('Reading Data')
data_g = xr.open_dataset(scr+'geopotential_500.nc')

#pdb.set_trace()
###So lets try to work this out
try:
    heights_smoothed = xr.open_dataset(mdir+'TOP_net_LW_all_daymean.nc')
except: 
    heights = data_g['z'].rename('heights')
    # Smooth the source data and re-insert into a new xarray DataArray
    heights_smoothed = xr.DataArray(savgol_filter(heights, 39, 2, mode="wrap"),coords=heights.coords,attrs=heights.attrs,name='heights_smoothed')
        
##################### Wave properties
# Variables for storage of analyzed values (lists of lists)
minima_full = [[] for lat in heights_smoothed.latitude]
maxima_full = [[] for lat in heights_smoothed.latitude]
wave_number_full = [[] for lat in heights_smoothed.latitude]
wave_amplitude_full = [[] for lat in heights_smoothed.latitude]
wave_speed_full = [[] for lat in heights_smoothed.latitude]

print('Starting Wave calculations! takes time!')
# Loop over hemispheres and times to compute wave properties (extrema, wave number, wave amplitude, wave speed)
for j in range(0, len(heights_smoothed.latitude)):
    print('Latitude: '+str(j)+' from '+str(len(heights_smoothed.latitude)))
    lat = heights_smoothed.latitude[j]
    # Get this particular latitude's storage variables
    minima = minima_full[j]
    maxima = maxima_full[j]
    wave_number = wave_number_full[j]
    wave_amplitude = wave_amplitude_full[j]
    wave_speed = wave_speed_full[j]
    for i in range(0, len(heights_smoothed.time)):
        # Filter out this time's height data
        current_heights = heights_smoothed.sel(latitude=lat).isel(time=i)
        current_heights = current_heights.heights_smoothed       
        # Get the minima and maxima
        minima_indicies = get_minima_indicies(current_heights)
        maxima_indicies = get_maxima_indicies(current_heights)
        minima.append(current_heights.longitude[minima_indicies].values)
        maxima.append(current_heights.longitude[maxima_indicies].values)
        # Calculate Wave Number
        k = max(len(minima_indicies), len(maxima_indicies))
        wave_number.append(k)
        # Calculate Amplitude
        trough_avg = current_heights[minima_indicies].values.mean()
        ridge_avg = current_heights[maxima_indicies].values.mean()
        a = (ridge_avg-trough_avg)/2.
        wave_amplitude.append(a)
        # Calculate Wave Speed
        if i == 0:
            # First time, cannot compute one-sided backwards difference
            c = np.nan
        else:
            # Otherwise, calculate speed (in deg lon/day)
            previous_minima = minima[i-1]
            previous_maxima = maxima[i-1]
            current_minima = minima[i]
            current_maxima = maxima[i]
            time_diff = heights_smoothed.time[i] - heights_smoothed.time[i-1]
            minima_increments = get_increments(previous_minima, current_minima)
            maxima_increments = get_increments(previous_maxima, current_maxima)
            lon_increments = np.concatenate([minima_increments, maxima_increments])
            c = reject_outliers(lon_increments).mean() / fraction_of_day(time_diff.values)
        wave_speed.append(c)

pdb.set_trace()
# Get our variables back into nice DataArrays
wave_property_coords = [('lat', heights_smoothed.latitude), ('time', heights_smoothed.time)]
wave_number_raw = xr.DataArray(wave_number_full, coords=wave_property_coords, name='wave_number')
wave_amplitude = xr.DataArray(wave_amplitude_full, coords=wave_property_coords, name='wave_amplitude')
wave_speed = xr.DataArray(wave_speed_full, coords=wave_property_coords, name='wave_speed')

# Set unit metadata
wave_number_raw.attrs['units'] = 'earth_radius**-1'
wave_amplitude.attrs['units'] = 'm'
wave_speed.attrs['units'] = 'degrees_east day**-1'

# Data to netcdf! finally a file!
wave_number_raw.to_netcdf(scr+'Wave_number.nc')
wave_amplitude.to_netcdf(scr+'Wave_amplitude.nc')
wave_speed.to_netcdf(scr+'Wave_speed.nc')

