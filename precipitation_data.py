import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import scipy.integrate as integrate
import scipy.special as special
import scipy.stats as stats
import os

# Files for MH:

file_aprl_MH= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/aprl_area_cut.nc"
data_aprl_MH= xr.open_dataset(file_aprl_MH)

file_aprc_MH= "C:/cygwin64/home/Ardhra/cdo/cut_area/MH/aprc_area_cut.nc"
data_aprc_MH= xr.open_dataset(file_aprc_MH)

# Files for PI:

file_aprl_PI= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/aprl_area_cut.nc"
data_aprl_PI= xr.open_dataset(file_aprl_PI)

file_aprc_PI= "C:/cygwin64/home/Ardhra/cdo/cut_area/PI/aprc_area_cut.nc"
data_aprc_PI= xr.open_dataset(file_aprc_PI)

file_name="//134.2.5.43/esd01/data/asmadhavan/scratch/Main_images/fig_"
my_path=os.path.abspath(file_name)

# Extracting variables for PI:

c_prec_PI= data_aprc_PI['aprc'][:]*60*60*24*30 # Converting to mm/month
l_prec_PI= data_aprl_PI['aprl'][:]*60*60*24*30


# Extracting variables for MH: 

c_prec_MH= data_aprc_MH['aprc'][:]*60*60*24*30
l_prec_MH= data_aprl_MH['aprl'][:]*60*60*24*30

# Separating seasons (Austral Summer and winter) PI:

c_prec_PI_JJA= c_prec_PI[5:8]                        # Take out the data for the months June, July, August (Austral winter)
c_prec_PI_JJA_mean= c_prec_PI_JJA.mean(dim='time')

l_prec_PI_JJA= l_prec_PI[5:8]
l_prec_PI_JJA_mean= l_prec_PI_JJA.mean(dim='time')

c_prec_PI_DJF= xr.concat([c_prec_PI[11],c_prec_PI[0],c_prec_PI[1]], dim='time') # Take out the data for December, Jan, and Feb and then concatenate them (Austral summer)
c_prec_PI_DJF_mean= c_prec_PI_DJF.mean(dim='time')

l_prec_PI_DJF= xr.concat([l_prec_PI[11],l_prec_PI[0],l_prec_PI[1]], dim='time')
l_prec_PI_DJF_mean= l_prec_PI_DJF.mean(dim='time')

prec_PI_JJA= c_prec_PI_JJA + l_prec_PI_JJA
prec_PI_JJA_mean= c_prec_PI_JJA_mean + l_prec_PI_JJA_mean   # Total precipitation is the sum of large-scale and convective precipitation

prec_PI_DJF= l_prec_PI_DJF + c_prec_PI_DJF
prec_PI_DJF_mean= l_prec_PI_DJF_mean + c_prec_PI_DJF_mean

# Separating seasons (Austral summer and winter) for MH:

c_prec_MH_JJA= c_prec_MH[5:8]
c_prec_MH_JJA_mean= c_prec_MH_JJA.mean(dim='time')

c_prec_MH_DJF= xr.concat([c_prec_MH[11], c_prec_MH[0], c_prec_MH[1]], dim='time')
c_prec_MH_DJF_mean= c_prec_MH_DJF.mean(dim='time')

l_prec_MH_JJA= l_prec_MH[5:8]
l_prec_MH_JJA_mean= l_prec_MH_JJA.mean(dim='time')

l_prec_MH_DJF= xr.concat([l_prec_MH[11],l_prec_MH[0],l_prec_MH[1]], dim='time')
l_prec_MH_DJF_mean= l_prec_MH_DJF.mean(dim='time')

prec_MH_JJA= c_prec_MH_JJA + l_prec_MH_JJA
prec_MH_JJA_mean= c_prec_MH_JJA_mean + l_prec_MH_JJA_mean

prec_MH_DJF= c_prec_MH_DJF + l_prec_MH_DJF
prec_MH_DJF_mean= c_prec_MH_DJF_mean + l_prec_MH_DJF_mean

diff_JJA= prec_MH_JJA_mean - prec_PI_JJA_mean
diff_DJF= prec_MH_DJF_mean - prec_PI_DJF_mean

# Statistical tests:

# Welch       : Gives a few significant values. I guess it gives more significant values for DJF.
# Ranksum test: Gives a few significant values

# Rank-sums p-values for JJA:

statistics_JJA, p_value_JJA= stats.ranksums(prec_PI_JJA,prec_MH_JJA) 

limit= 0.05                            # Setting the threshold value
mask= p_value_JJA < limit
pvalue_new_JJA= np.where(mask,p_value_JJA, np.nan)

# Welch's p-values for JJA:

statistics_JJA, p_value_JJA= stats.ttest_ind(prec_PI_JJA,prec_MH_JJA, equal_var=False) 

limit= 0.05                            # Setting the threshold value
mask= p_value_JJA < limit
pvalue_new_JJA= np.where(mask,p_value_JJA, np.nan)

# Rank-sums p-values for DJF:

statistics_DJF, p_value_DJF= stats.ranksums(prec_PI_DJF,prec_MH_DJF) 
limit= 0.05
mask= p_value_DJF < limit
pvalue_new_DJF= np.where(mask,p_value_DJF, np.nan)

# Welch's p-values for DJF:

statistics_DJF, p_value_DJF= stats.ttest_ind(prec_PI_DJF,prec_MH_DJF, equal_var=False) 
limit= 0.05
mask= p_value_DJF < limit
pvalue_new_DJF= np.where(mask,p_value_DJF, np.nan)


# Plotting:

fig, (ax1,ax2)= plt.subplots(1,2,figsize=(14,7),  subplot_kw={'projection':ccrs.PlateCarree()}) 
extent = [-80, -50, 0, -34 ]

lons= c_prec_MH['lon'][:]
lats= c_prec_MH['lat'][:]

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

plot1=ax1.pcolormesh(X,Y, prec_PI_JJA_mean, cmap='YlGnBu', transform=ccrs.PlateCarree(), vmin=0, vmax=600)
ax1.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax1.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax1.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot2=ax2.pcolormesh(X,Y, prec_PI_DJF_mean, cmap='YlGnBu', transform=ccrs.PlateCarree(), vmin=0, vmax=600)
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
cbar1.set_label(label='Precipitation [mm/month]', size= 25) # To change the size of the label [add weight='bold' for bold]
cbar1.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "AsEGUprec_ano_PI.svg", dpi=300)



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

plot3=ax3.pcolormesh(X,Y, diff_JJA, cmap='BrBG', transform=ccrs.PlateCarree(), vmin=-80, vmax=70)
ax3.scatter(X,Y, pvalue_new_JJA**0.5,facecolors='black', edgecolor='black', transform=ccrs.PlateCarree())
ax3.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
ax3.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
ax3.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
plot4=ax4.pcolormesh(X,Y, diff_DJF, cmap='BrBG', transform=ccrs.PlateCarree(), vmin=-80, vmax=70)
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


cbar2=fig.colorbar(plot4, shrink=1, ax=ax4, extend='both')
cbar2.set_label(label='Precipitation anomalies [mm/month]', size= 25) # To change the size of the label [add weight='bold']
cbar2.ax.tick_params(labelsize=22)  # Can be used to change the size of the values in the colorbar

plt.savefig(my_path + "AsEGUprec_ano_MH-PI_ranksums_withoutstat.svg", dpi=300)
plt.subplots_adjust(left=0.05, right=0.88, hspace=0.1)
plt.show()




