import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import plotly.figure_factory as ff
import scipy.stats as stats
import os


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

u_200= data_MH_u.sel(plev=20000)
v_200= data_MH_v.sel(plev=20000)
geopot_200= data_MH_geopot.sel(plev=20000)

lons= u_200['lon'][:]
lats= u_200['lat'][:]


# Data summer (DJF):

u_200_D= u_200[11]
u_200_J= u_200[0]
u_200_F= u_200[1]

u_200_DJF= xr.concat([u_200_D,u_200_J,u_200_F], dim='time').mean(dim='time')

# Data winter:
u200_JJA= u_200[5:8].mean(dim='time')



# Data summer V (DJF):

v_200_D= v_200[11]
v_200_J= v_200[0]
v_200_F= v_200[1]

v_200_DJF= xr.concat([v_200_D,v_200_J,v_200_F], dim='time').mean(dim='time')


# Data winter:
v200_JJA= v_200[5:8].mean(dim='time')




# Create coordinate grid

# print(u_200_DJF)

#############################

# Data summer (DJF):

geopot_200_D= geopot_200[11]
geopot_200_J= geopot_200[0]
geopot_200_F= geopot_200[1]

geopot_200_DJF= xr.concat([geopot_200_D,geopot_200_J,geopot_200_F], dim='time').mean(dim='time')

# Data winter:
geopot_200_JJA= geopot_200[5:8].mean(dim='time')

###################################

########### DATA PI ########################
data_PI_u= data_PI_u['u']
data_PI_v=data_PI_v['v']
data_PI_geopot= data_PI_geopot['geopoth']

u_200_PI= data_PI_u.sel(plev=20000)
v_200_PI= data_PI_v.sel(plev=20000)
geopot_200_PI= data_PI_geopot.sel(plev=20000)

lons= u_200_PI['lon'][:]
lats= u_200_PI['lat'][:]


# Data summer (DJF):

u_200_PI_D= u_200_PI[11]
u_200_PI_J= u_200_PI[0]
u_200_PI_F= u_200_PI[1]

u_200_PI_DJF= xr.concat([u_200_PI_D,u_200_PI_J,u_200_PI_F], dim='time').mean(dim='time')

# Data winter:
u200_PI_JJA= u_200_PI[5:8].mean(dim='time')



# Data summer V (DJF):

v_200_PI_D= v_200_PI[11]
v_200_PI_J= v_200_PI[0]
v_200_PI_F= v_200_PI[1]

v_200_PI_DJF= xr.concat([v_200_PI_D,v_200_PI_J,v_200_PI_F], dim='time').mean(dim='time')


# Data winter:
v200_PI_JJA= v_200_PI[5:8].mean(dim='time')

# Data geopot DJF

#print(geopot_200)

geopot_200_PI_D= geopot_200_PI[11]
geopot_200_PI_J= geopot_200_PI[0]
geopot_200_PI_F= geopot_200_PI[1]

geopot_200_PI_DJF= xr.concat([geopot_200_PI_D,geopot_200_PI_J,geopot_200_PI_F], dim='time').mean(dim='time')

# Data winter:
geopot_200_PI_JJA= geopot_200_PI[5:8].mean(dim='time')

# Create coordinate grid

# print(u_200_DJF)


###########################################################################

u1_djf= u_200_DJF - u_200_PI_DJF
v1_djf= v_200_DJF -v_200_PI_DJF

u1_jja= u200_JJA - u200_PI_JJA
v1_jja= v200_JJA - v200_PI_JJA

geopot_ano_djf= geopot_200_DJF - geopot_200_PI_DJF
geopot_ano_jja= geopot_200_JJA - geopot_200_PI_JJA


X=lons
Y=lats
X, Y= np.meshgrid(X,Y)

#################################################################################################################################################################################################
##########################################################################################################################

# Files MH:

file_aprc_MH= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/aprc_area_cut.nc"
data_aprc_MH= xr.open_dataset(file_aprc_MH)

file_aprl_MH= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/aprl_area_cut.nc"
data_aprl_MH= xr.open_dataset(file_aprl_MH)

# Files PI:

file_aprc_PI="C:/cygwin64/home/Ardhra/cdo/cut_area/PI/aprc_area_cut.nc"
data_aprc_PI= xr.open_dataset(file_aprc_PI)

file_aprl_PI= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/aprl_area_cut.nc"
data_aprl_PI= xr.open_dataset(file_aprl_PI)

# aprl DJF MH:

aprl_MH= data_aprl_MH['aprl'][:]
lon= aprl_MH['lon'][:]
lat= aprl_MH['lat'][:]
aprl_MH_DJF= xr.concat([aprl_MH[11], aprl_MH[0], aprl_MH[1]], dim='time')

# aprl DJF PI: 

aprl_PI= data_aprl_PI['aprl'][:]

aprl_PI_DJF= xr.concat([aprl_PI[11], aprl_PI[0],aprl_PI[1]], dim='time')

# aprc DJF MH:

aprc_MH= data_aprc_MH['aprc'][:]

aprc_MH_DJF= xr.concat([aprc_MH[11], aprc_MH[0], aprc_MH[1]], dim='time')

# aprc DJF PI:

aprc_PI= data_aprc_PI['aprc'][:]

aprc_PI_DJF= xr.concat([aprc_PI[11],aprc_PI[0],aprc_PI[1]], dim='time')

# Prec DJF PI:

prec_PI_DJF= aprc_PI_DJF + aprl_PI_DJF

# Prec DJF MH:

prec_MH_DJF= aprc_MH_DJF + aprl_MH_DJF

# Datasets

data_1= prec_PI_DJF*60*60*24*30

data_2= (prec_MH_DJF - prec_PI_DJF)*60*60*24*30

# Welch's test:

statistic, p_values = stats.ttest_ind(data_1, data_2, axis=0, equal_var=False)
#print(f"pvale:{np.shape(p_values)}")

limit= 0.01
# # Select values below the limit

mask= p_values < limit

pvalue_new = np.where(mask,p_values, np.nan)


#########################################################################################################################################################################################
##########################################################################################################################

# Create streamline plot
# fig, axs = plt.subplots(2,2,figsize=(15,13), subplot_kw={'projection':ccrs.PlateCarree()})

# ax1= axs[0][0]
# ax2= axs[0][1]
# ax3= axs[1][0]
# ax4= axs[1][1]

# ax1.add_feature(cfeature.LAND)
# ax1.add_feature(cfeature.COASTLINE, color='k')
# ax2.add_feature(cfeature.LAND)
# ax2.add_feature(cfeature.COASTLINE, color='k')
# ax1.set_title('PI (JJA)', loc='left', fontsize=25)
# ax2.set_title('PI (DJF)', loc='left', fontsize=25)
# ax3.add_feature(cfeature.LAND)
# ax3.add_feature(cfeature.COASTLINE, color='k')
# ax4.add_feature(cfeature.LAND)
# ax4.add_feature(cfeature.COASTLINE, color='k')
# ax3.set_title('MH-PI (JJA)', loc='left', fontsize=25)
# ax4.set_title('MH-PI (DJF)', loc='left', fontsize=25)

# #  Varying density along a streamline
# plot1= ax1.pcolormesh(X,Y, geopot_200_PI_JJA,cmap= 'RdBu', transform=ccrs.PlateCarree())
# ax1.quiver(X,Y, u200_PI_JJA, v200_PI_JJA, transform= ccrs.PlateCarree() , scale=160)
# plot2= ax2.pcolormesh(X,Y, geopot_200_PI_DJF,cmap= 'RdBu', transform=ccrs.PlateCarree())
# ax2.quiver(X,Y, u_200_PI_DJF, v_200_PI_DJF, transform= ccrs.PlateCarree() , scale=65)
# # ax2.scatter(X,Y, pvalue_new**0.5,facecolors='white', edgecolor='white', transform=ccrs.PlateCarree())
# plot3= ax3.pcolormesh(X,Y, geopot_200_JJA,cmap= 'RdBu', transform=ccrs.PlateCarree())
# ax3.quiver(X,Y, u200_JJA, v200_JJA, transform= ccrs.PlateCarree(), scale=160)
# plot4= ax4.pcolormesh(X,Y, geopot_200_DJF,cmap= 'RdBu', transform=ccrs.PlateCarree())
# ax4.quiver(X,Y, u_200_DJF, v_200_DJF, transform= ccrs.PlateCarree() , scale=65)
# # ax4.scatter(X,Y, pvalue_new**0.5,facecolors='white', edgecolor='white', transform=ccrs.PlateCarree())
# fig.tight_layout(pad=5.0)

# cbar1=fig.colorbar(plot2, shrink=1, ax=ax2, format='%.0f', extend='both')
# cbar1.set_label(label='Geopotential height [m]', size= 25) # To change the size of the label [add weight='bold' for bold italice whatever]
# cbar1.ax.tick_params(labelsize=25)  # Can be used to change the size of the values in the colorbar
# cbar2=fig.colorbar(plot4, shrink=1, ax=ax4, format='%.0f')
# cbar2.set_label(label='Geopotential height difference [m]', size= 25) # To change the size of the label [add weight='bold' for bold italice whatever]
# cbar2.ax.tick_params(labelsize=25)  # Can be used to change the size of the values in the colorbar


# fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
# plt.tight_layout()
# plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.10, hspace=0.1)
# fig.show()
# fig.savefig('//134.2.5.43/esd01/data/asmadhavan/scratch/Main_images/streamline_plot.svg', dpi=300, bbox_inches='tight')    

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


# plt.colorbar(plot1, shrink=1, orientation='vertical', label= 'Precipitation [mm/month]')
#cax = plt.axes([0.90, 0.11, 0.025, 0.775])
cbar1=fig.colorbar(plot2, shrink=1, ax=ax2, extend='both') #, format='%.0f'
cbar1.set_label(label='Geopotential height [m]', size= 25) # To change the size of the label [add weight='bold' for bold italice whatever]
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
cbar2.set_label(label='Geopotential height [m]', size= 25) # To change the size of the label [add weight='bold' for bold italice whatever]
cbar2.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "geopot_streamline_MH.svg", dpi=300)
plt.subplots_adjust(left=0.05, right=0.88, hspace=0.1)
plt.show()






