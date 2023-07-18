import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.stats as stats

# file MH: 

file_MH= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/rhumidity_levels_areacut.nc"
data_MH= xr.open_dataset(file_MH)

# File PI:

file_PI= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/rhumidity_levels_areacut.nc"
data_PI= xr.open_dataset(file_PI)

file_name="//134.2.5.43/esd01/data/asmadhavan/scratch/Main_images/fig_"
my_path=os.path.abspath(file_name)

# Extracting the variables:

relhum_MH= data_MH['relhum'][:] #[time, plevel, lon, lat]
relhum_PI= data_PI['relhum'][:]

minlev= 100000                                          # Selecting the minimum and maximum vertical levels to be selected.
maxlev= 20000

relhum_MH= relhum_MH.sel(plev= slice(minlev, maxlev))   # Slicing the data within the min and max vertical levels
relhum_PI= relhum_PI.sel(plev= slice(minlev, maxlev))

relhum_MH= relhum_MH.mean(dim='plev')                   # Averaging along the different pressure levels to find mean value.
relhum_PI= relhum_PI.mean(dim='plev')

lons= relhum_PI['lon'][:]
lats= relhum_PI['lat'][:]

# MH season separation:
relhum_MH_DJF= xr.concat([relhum_MH[11], relhum_MH[0], relhum_MH[1]], dim='time')                             # Austral summer
relhum_MH_DJF_plot= xr.concat([relhum_MH[11], relhum_MH[0], relhum_MH[1]], dim='time').mean(dim='time')
relhum_MH_JJA= relhum_MH[5:8]                                                                                 # Austral winter
relhum_MH_JJA_plot= relhum_MH[5:8].mean(dim='time')

# PI season separation:

relhum_PI_DJF= xr.concat([relhum_PI[11], relhum_PI[0], relhum_PI[1]], dim='time')
relhum_PI_DJF_plot= xr.concat([relhum_PI[11], relhum_PI[0], relhum_PI[1]], dim='time').mean(dim='time')
relhum_PI_JJA= relhum_PI[5:8]
relhum_PI_JJA_plot= relhum_PI[5:8].mean(dim='time')

# Differences (MH-PI):

relhum_MH_PI_DJF= relhum_MH_DJF_plot - relhum_PI_DJF_plot
relhum_MH_PI_JJA= relhum_MH_JJA_plot - relhum_PI_JJA_plot

# Rank-sums Statistical analysis JJA:

data_1_JJA= relhum_PI_JJA
data_2_JJA= relhum_MH_JJA
statistic_JJA, p_values_JJA = stats.ranksums(data_1_JJA, data_2_JJA)
limit=0.05
mask= p_values_JJA < limit
pvalue_new_JJA = np.where(mask,p_values_JJA, np.nan)

# Welch's statistical analysis JJA:
statistic_JJA, p_values_JJA = stats.ttest_ind(data_1_JJA, data_2_JJA, equal_var=False)
limit=0.05
mask= p_values_JJA < limit
pvalue_new_JJA = np.where(mask,p_values_JJA, np.nan)

# Rank-sums statistical analysis DJF:

data_1_DJF= relhum_PI_DJF
data_2_DJF= relhum_MH_DJF
statistic_DJF, p_values_DJF= stats.ranksums(data_1_DJF, data_2_DJF)
limit=0.05
mask= p_values_DJF < limit
pvalue_new_DJF= np.where(mask, p_values_DJF, np.nan)

# Welchs's statistical analysis DJF: 
statistic_DJF, p_values_DJF= stats.ttest_ind(data_1_DJF, data_2_DJF, equal_var=False)
limit=0.05
mask= p_values_DJF < limit
pvalue_new_DJF= np.where(mask, p_values_DJF, np.nan)

# Plotting:

fig, (ax1,ax2)= plt.subplots(1,2,figsize=(14,7),  subplot_kw={'projection':ccrs.PlateCarree()}) 
extent = [-80, -50, 0, -34 ]

X=lons
Y=lats
X,Y=np.meshgrid(X,Y)

ax1.add_feature(cfeature.LAND)
ax1.add_feature(cfeature.COASTLINE, color='k')
ax2.add_feature(cfeature.LAND)
ax2.add_feature(cfeature.COASTLINE, color='k')
ax1.set_extent(extent)
ax2.set_extent(extent)
ax1.set_title('PI (JJA)', loc='left', fontsize=25)
ax2.set_title('PI (DJF)', loc='left', fontsize=25)

plot1=ax1.pcolormesh(X,Y, relhum_PI_JJA_plot, cmap='RdBu', transform=ccrs.PlateCarree())
ax1.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax1.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot2=ax2.pcolormesh(X,Y, relhum_PI_DJF_plot, cmap='RdBu', transform=ccrs.PlateCarree())
ax2.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax2.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax2.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax2.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
fig.tight_layout(pad=5.0)

gl= ax1.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 25}
gl.ylabel_style = {'size': 25}


gl2= ax2.gridlines(draw_labels=True)
gl2.top_labels = False
gl2.right_labels = False
gl2.xlabel_style = {'size': 25}
gl2.ylabel_style = {'size': 25}

cbar1=fig.colorbar(plot2, shrink=1, ax=ax2, extend='neither')  #, format='%.0f'
cbar1.set_label(label='Relative humidity', size= 25) # To change the size of the label [add weight='bold']
cbar1.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "AsEGU_relhum_ranksums.svg", dpi=300)

fig, (ax3,ax4)= plt.subplots(1,2, figsize=(14,7), subplot_kw={'projection':ccrs.PlateCarree()})
extent = [-80, -50, 0, -34 ]
level= np.linspace(-100,80,35)

ax3.set_title('MH-PI (JJA)', loc='left', fontsize=25)
ax4.set_title('MH-PI (DJF)', loc='left', fontsize=25)
ax3.add_feature(cfeature.LAND)
ax3.add_feature(cfeature.COASTLINE, color='k')
ax4.add_feature(cfeature.LAND)
ax4.add_feature(cfeature.COASTLINE, color='k')
ax3.set_extent(extent)
ax4.set_extent(extent)

vmin=-0.04
vmax= 0.04

plot3=ax3.pcolormesh(X,Y, relhum_MH_PI_JJA, cmap='RdBu', transform=ccrs.PlateCarree(), vmin=-0.04, vmax= 0.04)
ax3.scatter(X,Y, pvalue_new_JJA**0.5,facecolors='black', edgecolor='black', transform=ccrs.PlateCarree())
ax3.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax3.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot4=ax4.pcolormesh(X,Y, relhum_MH_PI_DJF, cmap='RdBu', transform=ccrs.PlateCarree(), vmin=-0.04, vmax= 0.04)
ax4.scatter(X,Y, pvalue_new_DJF**0.5,facecolors='black', edgecolor='black', transform=ccrs.PlateCarree())
ax4.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax4.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax4.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax4.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
fig.tight_layout(pad=5.0)

gl3= ax3.gridlines(draw_labels=True)
gl3.top_labels = False
gl3.right_labels = False
gl3.xlabel_style = {'size': 25}
gl3.ylabel_style = {'size': 25}


gl4= ax4.gridlines(draw_labels=True)
gl4.top_labels = False
gl4.right_labels = False
gl4.xlabel_style = {'size': 25}
gl4.ylabel_style = {'size': 25}

cbar2=fig.colorbar(plot4, shrink=1, ax=ax4, extend='both') # , format='%.0f'
cbar2.set_label(label='Relative humidity', size= 25) # To change the size of the label [add weight='bold' for bold]
cbar2.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar


plt.savefig(my_path + "AsEGU_MH-PI_relhum_ranksums.svg", dpi=300)
plt.subplots_adjust(left=0.05, right=0.88, hspace=0.1)
plt.show()


