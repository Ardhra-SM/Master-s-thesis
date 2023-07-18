import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import scipy.integrate as integrate
import scipy.special as special
import scipy.stats as stats
import os

file_PI_q= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/specific_humidity_levels_areacut.nc"
data_PI_q= xr.open_dataset(file_PI_q)

file_MH_q= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/specific_humidity_levels_areacut.nc"
data_MH_q= xr.open_dataset(file_MH_q)

file_PI_u= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/u_levels_areacut.nc"
data_PI_u= xr.open_dataset(file_PI_u)

file_MH_u= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/u_levels_areacut.nc"
data_MH_u= xr.open_dataset(file_MH_u)

file_PI_v= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/v_levels_areacut.nc"
data_PI_v= xr.open_dataset(file_PI_v)

file_MH_v= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/v_levels_areacut.nc"
data_MH_v= xr.open_dataset(file_MH_v)

file_MH_geopot= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/geopot_area_cut.nc"
data_MH_geopot= xr.open_dataset(file_MH_geopot)

file_PI_geopot= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/geopoth_area_cut.nc"
data_PI_geopot= xr.open_dataset(file_PI_geopot)

file_name="//134.2.5.43/esd01/data/asmadhavan/scratch/Main_images/fig_"
my_path=os.path.abspath(file_name)

# Separating the variables:
min_plev= 100000
max_plev= 20000

lons= data_MH_u['lon'][:]
lats= data_MH_u['lat'][:]


u_PI= data_PI_u['u']                          # u-wind data
u_PI= u_PI.sel(plev=slice(min_plev,max_plev)) # To select the whole range of plev
u_MH= data_MH_u['u']
u_MH= u_MH.sel(plev=slice(min_plev,max_plev))

v_PI= data_PI_v['v']                          # v-win data              
v_PI=v_PI.sel(plev=slice(min_plev,max_plev))
v_MH= data_MH_v['v']
v_MH=v_MH.sel(plev=slice(min_plev,max_plev))

q_PI= data_PI_q['q']                          # Specific humidity data
q_PI= q_PI.sel(plev=slice(min_plev,max_plev))
q_MH= data_MH_q['q']
q_MH= q_MH.sel(plev=slice(min_plev,max_plev))

geopot_PI= data_PI_geopot['geopoth']          # Geopotential height data
geopot_PI= geopot_PI.sel(plev=slice(min_plev,max_plev))
geopot_MH= data_MH_geopot['geopoth']
geopot_MH= geopot_MH.sel(plev=slice(min_plev,max_plev))

# Data 'u' for PI:

u_PI_JJA= u_PI[5:8]
u_PI_JJA_plot= np.mean(u_PI_JJA, axis=0)

u_PI_DJF= xr.concat([u_PI[11],u_PI[0],u_PI[1]], dim='time')
u_PI_DJF_plot= u_PI_DJF.mean(dim='time')

# Data 'v' for PI:

v_PI_JJA= v_PI[5:8]
v_PI_JJA_plot= np.mean(v_PI_JJA, axis=0)

v_PI_DJF= xr.concat([v_PI[11],v_PI[0],v_PI[1]], dim='time')
v_PI_DJF_plot= v_PI_DJF.mean(dim='time')

# Data 'q' for PI:

q_PI_JJA= q_PI[5:8]
q_PI_JJA_plot= np.mean(q_PI_JJA, axis=0)

q_PI_DJF= xr.concat([q_PI[11], q_PI[0],q_PI[1]], dim='time')
q_PI_DJF_plot= q_PI_DJF.mean(dim='time')

# Data geopot for PI:

geopot_PI_JJA= geopot_PI[5:8]
geopot_PI_JJA_plot= np.mean(geopot_PI_JJA, axis=0)

geopot_PI_DJF= xr.concat([geopot_PI[11], geopot_PI[0],geopot_PI[1]], dim='time')
geopot_PI_DJF_plot= geopot_PI_DJF.mean(dim='time')

# Data 'u' for MH:

u_MH_JJA= u_MH[5:8]
u_MH_JJA_plot= u_MH_JJA.mean(dim='time')

u_MH_DJF= xr.concat([u_MH[11], u_MH[0],u_MH[1]], dim='time')
u_MH_DJF_plot= u_MH_DJF.mean(dim='time')

# Data 'v' for MH:

v_MH_JJA= v_MH[5:8]
v_MH_JJA_plot= v_MH_JJA.mean(dim='time')

v_MH_DJF= xr.concat([v_MH[11], v_MH[0], v_MH[1]], dim='time')
v_MH_DJF_plot= v_MH_DJF.mean(dim='time')

# Data 'q' for MH:

q_MH_JJA= q_MH[5:8]
q_MH_JJA_plot= q_MH_JJA.mean(dim='time')

q_MH_DJF= xr.concat([q_MH[11], q_MH[0],q_MH[1]], dim='time')
q_MH_DJF_plot= q_MH_DJF.mean(dim='time')

# Data geopot for MH:

geopot_MH_JJA= geopot_MH[5:8]
geopot_MH_JJA_plot=geopot_MH_JJA.mean(dim='time')

geopot_MH_DJF= xr.concat([geopot_MH[11], geopot_MH[0], geopot_MH[1]], dim='time')
geopot_MH_DJF_plot= geopot_MH_DJF.mean(dim='time')

################### VIMFC for statistical analysis #######################


vimfc_PI_JJA_stats= ((q_PI_JJA*u_PI_JJA)+(q_PI_JJA*v_PI_JJA))/-9.81
vimfc_PI_JJA_stats= vimfc_PI_JJA_stats.integrate(coord='plev')

vimfc_PI_DJF_stats= ((q_PI_DJF*u_PI_DJF)+(q_PI_DJF*v_PI_DJF))/-9.81
vimfc_PI_DJF_stats= vimfc_PI_DJF_stats.integrate(coord='plev')


vimfc_MH_JJA_stats= ((q_MH_JJA*u_MH_JJA)+(q_MH_JJA*v_MH_JJA))/-9.81
vimfc_MH_JJA_stats= vimfc_MH_JJA_stats.integrate(coord='plev')


vimfc_MH_DJF_stats= ((q_MH_DJF*u_MH_DJF)+(q_MH_DJF*v_MH_DJF))/-9.81
vimfc_MH_DJF_stats= vimfc_MH_DJF_stats.integrate(coord='plev')

################################## VIMF Zonal component ###################

vimfc_PI_JJA_zonal= (q_PI_JJA*u_PI_JJA)/-9.81
vimfc_PI_JJA_zonal= vimfc_PI_JJA_zonal.integrate(coord='plev')
vimfc_PI_JJA_zonal= vimfc_PI_JJA_zonal.mean(dim='time')

vimfc_PI_DJF_zonal= (q_PI_DJF*u_PI_DJF)/-9.81
vimfc_PI_DJF_zonal= vimfc_PI_DJF_zonal.integrate(coord='plev')
vimfc_PI_DJF_zonal= vimfc_PI_DJF_zonal.mean(dim='time')

vimfc_MH_JJA_zonal= (q_MH_JJA*u_MH_JJA)/-9.81
vimfc_MH_JJA_zonal= vimfc_MH_JJA_zonal.integrate(coord='plev')
vimfc_MH_JJA_zonal= vimfc_MH_JJA_zonal.mean(dim='time')

vimfc_MH_DJF_zonal= (q_MH_DJF*u_MH_DJF)/-9.81
vimfc_MH_DJF_zonal= vimfc_MH_DJF_zonal.integrate(coord='plev')
vimfc_MH_DJF_zonal= vimfc_MH_DJF_zonal.mean(dim='time')

################################## VIMF Meridional component ###################

vimfc_PI_JJA_meridonal= (q_PI_JJA*v_PI_JJA)/-9.81
vimfc_PI_JJA_meridonal= vimfc_PI_JJA_meridonal.integrate(coord='plev')
vimfc_PI_JJA_meridonal= vimfc_PI_JJA_meridonal.mean(dim='time')

vimfc_PI_DJF_meridional= (q_PI_DJF*v_PI_DJF)/-9.81
vimfc_PI_DJF_meridional= vimfc_PI_DJF_meridional.integrate(coord='plev')
vimfc_PI_DJF_meridional= vimfc_PI_DJF_meridional.mean(dim='time')

vimfc_MH_JJA_meridional= (q_MH_JJA*v_MH_JJA)/-9.81
vimfc_MH_JJA_meridional= vimfc_MH_JJA_meridional.integrate(coord='plev')
vimfc_MH_JJA_meridional= vimfc_MH_JJA_meridional.mean(dim='time')

vimfc_MH_DJF_meridional= (q_MH_DJF*v_MH_DJF)/-9.81
vimfc_MH_DJF_meridional= vimfc_MH_DJF_meridional.integrate(coord='plev')
vimfc_MH_DJF_meridional= vimfc_MH_DJF_meridional.mean(dim='time')

############################ VIMFC for plots ##############################

vimfc_PI_JJA= np.mean(vimfc_PI_JJA_stats, axis=0)
vimfc_PI_DJF= np.mean(vimfc_PI_DJF_stats, axis=0)
vimfc_MH_JJA= np.mean(vimfc_MH_JJA_stats, axis=0)
vimfc_MH_DJF= np.mean(vimfc_MH_DJF_stats, axis=0)

vimfc_MH_PI_JJA= vimfc_MH_JJA - vimfc_PI_JJA
vimfc_MH_PI_DJF= vimfc_MH_DJF - vimfc_PI_DJF

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

plot1=ax1.pcolormesh(X,Y, vimfc_PI_JJA, cmap='coolwarm', transform=ccrs.PlateCarree(), vmin=-250, vmax=250)
ax1.quiver(X[::2],Y[::2], vimfc_PI_JJA_zonal[::2], vimfc_PI_JJA_meridonal[::2], transform= ccrs.PlateCarree() , color='black', scale=1050) # , scale=350
ax1.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax1.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot2=ax2.pcolormesh(X,Y, vimfc_PI_DJF, cmap='coolwarm', transform=ccrs.PlateCarree(), vmin=-250, vmax=250)
ax2.quiver(X[::2],Y[::2], vimfc_PI_DJF_zonal[::2], vimfc_PI_DJF_meridional[::2], transform= ccrs.PlateCarree() , color='black', scale=1050)
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
cbar1.set_label(label='VIMF [kg/ms]', size= 25) # To change the size of the label [add weight='bold' for bold]
cbar1.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "AsEGUvimfc_ano_PI_ranksums_0.1thresh.svg", dpi=300)

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

plot3=ax3.pcolormesh(X,Y, vimfc_MH_PI_JJA , cmap='coolwarm', transform=ccrs.PlateCarree(), vmin=-30, vmax=30)
ax3.quiver(X[::2],Y[::2], vimfc_MH_JJA_zonal[::2] - vimfc_PI_JJA_zonal[::2], vimfc_MH_JJA_meridional[::2] - vimfc_PI_JJA_meridonal[::2], transform= ccrs.PlateCarree() , color='black', scale=350)
ax3.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax3.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot4=ax4.pcolormesh(X,Y, vimfc_MH_PI_DJF, cmap='coolwarm', transform=ccrs.PlateCarree(), vmin=-30, vmax=30)
ax4.quiver(X[::2],Y[::2], vimfc_MH_DJF_zonal[::2] - vimfc_PI_DJF_zonal[::2], vimfc_MH_DJF_meridional[::2] - vimfc_PI_DJF_meridional[::2], transform= ccrs.PlateCarree() , color='black', scale=350)
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
cbar2.set_label(label='VIMF [kg/ms]', size= 25) # To change the size of the label [add weight='bold' for bold]
cbar2.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "AsEGUvimfc_ano_MH-PI_ranksums_200623_quiver.svg", dpi=300)
plt.subplots_adjust(left=0.05, right=0.88, hspace=0.1)
plt.show()
