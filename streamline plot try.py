import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import plotly.figure_factory as ff
import scipy.stats as stats
import os

# File MH:

file_MH_u= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/u_levels_areacut.nc"
data_MH_u= xr.open_dataset(file_MH_u)

file_MH_v= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/v_levels_areacut.nc"
data_MH_v= xr.open_dataset(file_MH_v)

file_MH_geopot= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/geopot_area_cut.nc"
data_MH_geopot= xr.open_dataset(file_MH_geopot)

# File PI:

file_PI_u= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/u_levels_areacut.nc"
data_PI_u= xr.open_dataset(file_PI_u)

file_PI_v= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/v_levels_areacut.nc"
data_PI_v= xr.open_dataset(file_PI_v)

file_PI_geopot= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/geopoth_area_cut.nc"
data_PI_geopot= xr.open_dataset(file_PI_geopot)

file_name="//134.2.5.43/esd01/data/asmadhavan/scratch/Main_images/fig_"
my_path=os.path.abspath(file_name)

data_MH_u= data_MH_u['u']
data_MH_v=data_MH_v['v']
data_MH_geopot= data_MH_geopot['geopoth']

# Selecting the vertical level for each of the variables
u_200= data_MH_u.sel(plev=20000)
v_200= data_MH_v.sel(plev=20000)
geopot_200= data_MH_geopot.sel(plev=20000)

lons= u_200['lon'][:]
lats= u_200['lat'][:]


# MH season separation for u-wind:

u_200_DJF= xr.concat([u_200[11],u_200[0],u_200[1]], dim='time').mean(dim='time')  # Austral summer
u200_JJA= u_200[5:8].mean(dim='time')                                             # Austral winter

# MH season separation for v-wind:

v_200_DJF= xr.concat([v_200[11],v_200[0],v_200[1]], dim='time').mean(dim='time')
v200_JJA= v_200[5:8].mean(dim='time')

# MH season separation for geopotential height:

geopot_200_DJF= xr.concat([geopot_200[11], geopot_200[0],geopot_200[1]], dim='time').mean(dim='time')
geopot_200_JJA= geopot_200[5:8].mean(dim='time')


# Data PI
data_PI_u= data_PI_u['u']
data_PI_v=data_PI_v['v']
data_PI_geopot= data_PI_geopot['geopoth']

u_200_PI= data_PI_u.sel(plev=20000)
v_200_PI= data_PI_v.sel(plev=20000)
geopot_200_PI= data_PI_geopot.sel(plev=20000)

lons= u_200_PI['lon'][:]
lats= u_200_PI['lat'][:]


# PI season separation for u-wind:

u_200_PI_DJF= xr.concat([u_200_PI[11],u_200_PI[0],u_200_PI[1]], dim='time').mean(dim='time')
u200_PI_JJA= u_200_PI[5:8].mean(dim='time')

# PI season separation for v-wind:

v_200_PI_DJF= xr.concat([v_200_PI[11],v_200_PI[0],v_200_PI[1]], dim='time').mean(dim='time')
v200_PI_JJA= v_200_PI[5:8].mean(dim='time')

# PI season separation of geopotential height:

geopot_200_PI_DJF= xr.concat([geopot_200_PI[11],geopot_200_PI[0],geopot_200_PI[1]], dim='time').mean(dim='time')
geopot_200_PI_JJA= geopot_200_PI[5:8].mean(dim='time')

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

plot1=ax1.pcolormesh(X,Y, geopot_200_PI_JJA, cmap='RdBu', transform=ccrs.PlateCarree())
ax1.quiver(X[::2],Y[::2], u200_PI_JJA[::2], v200_PI_JJA[::2], transform= ccrs.PlateCarree() , color='white', scale=350)
ax1.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax1.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot2=ax2.pcolormesh(X,Y, geopot_200_PI_DJF, cmap='RdBu', transform=ccrs.PlateCarree())
ax2.quiver(X[::2],Y[::2], u_200_PI_DJF[::2], v_200_PI_DJF[::2], transform= ccrs.PlateCarree() , color='white' , scale=100)
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

cbar1=fig.colorbar(plot2, shrink=1, ax=ax2, extend='both') #, format='%.0f'
cbar1.set_label(label='Geopotential height [m]', size= 25) # To change the size of the label [add weight='bold' for bold]
cbar1.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "geopot_streamline_PI.svg", dpi=300)

fig, (ax3,ax4)= plt.subplots(1,2, figsize=(14,7), subplot_kw={'projection':ccrs.PlateCarree()})
extent = [-80, -50, 0, -34 ]
level= np.linspace(-100,80,35)

ax3.set_title('MH (JJA)', loc='left', fontsize=25)
ax4.set_title('MH (DJF)', loc='left', fontsize=25)
ax3.add_feature(cfeature.LAND)
ax3.add_feature(cfeature.COASTLINE, color='k')
ax4.add_feature(cfeature.LAND)
ax4.add_feature(cfeature.COASTLINE, color='k')
ax3.set_extent(extent)
ax4.set_extent(extent)

plot3=ax3.pcolormesh(X,Y, geopot_200_JJA, cmap='RdBu', transform=ccrs.PlateCarree())
ax3.quiver(X[::2],Y[::2], u200_JJA[::2], v200_JJA[::2], transform= ccrs.PlateCarree(), color= 'white', scale=350)
ax3.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax3.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot4=ax4.pcolormesh(X,Y, geopot_200_DJF, cmap='RdBu', transform=ccrs.PlateCarree())
ax4.quiver(X[::2],Y[::2], u_200_DJF[::2], v_200_DJF[::2], transform= ccrs.PlateCarree() , color='white' ,scale=100)
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


cbar2=fig.colorbar(plot4, shrink=1, ax=ax4, extend='both')
cbar2.set_label(label='Geopotential height [m]', size= 25) # To change the size of the label [add weight='bold' for bold]
cbar2.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "geopot_streamline_MH.svg", dpi=300)
plt.subplots_adjust(left=0.05, right=0.88, hspace=0.1)
plt.show()






