# have a dataset
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

# Converts precipitation unit from kg/m2 to mm/month
def prec_mm_month(prec):
    mm_month= prec*60*60*24*30
    print("----------Precipitation has been converted to mm/month----------")
    return mm_month
--------------------------------------------------------------------------------------------------------------------------
# Extracts the "convective precipitation" variable from the netCDF file
def convective_precipitation(file_name):
    c_prec= file_name['aprc'][:]
    print("----------convective precipitation has been extracted----------")
    return c_prec
--------------------------------------------------------------------------------------------------------------------------------
# Extracts the "large scale precipitation" variable from the netCDF file 
def largescale_precipitation(file_name):
    l_prec= file_name['aprl'][:]
    print("----------large scale precipitation has been extracted----------")
    return l_prec
--------------------------------------------------------------------------------------------------------------------------------
# Extracts only the austral summer months
def summer_south(data):
    sum_south= xr.concat([data[11], data[0], data[1]], dim='time')
    print("----------Summer south data created----------")
    return sum_south
-----------------------------------------------------------------------------------------------------------------------------
# Extracts only the austral winter months
def winter_south(data):
    win_south= data[5:8]
    print("----------Winter south data created----------")
    return win_south
-------------------------------------------------------------------------------------------------------------------------------
# Extracts the total precipitation from convective and large-scale precipitation
def total_prec(convective, large):
    total= convective + large
    print("----------Total precipitation has been calculated----------")
    return total
-------------------------------------------------------------------------------------------------------------------------------
# Extracts the longitude values from the data
def longitude(filename):
    lon= filename['lon'][:]
    return lon
------------------------------------------------------------------------------------------------------------    
# Extracts the latitude values from the data
def latitude(filename):
    lat= filename['lat'][:]
    return lat
--------------------------------------------------------------------------------------------------------------
# Plots a 2 subplot figure when provided with var1, var2, lon, lat, and the extent
def two_fig(varname1, varname2, lons, lats, extent):
    fig, (ax1,ax2)= plt.subplots(1,2,figsize=(14,7),  subplot_kw={'projection':ccrs.PlateCarree()}) 
    # extent = [-80, -50, 0, -34 ]
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

    plot1=ax1.pcolormesh(X,Y, varname1, cmap='YlGnBu', transform=ccrs.PlateCarree(), vmin=0, vmax=600)
    plot2=ax2.pcolormesh(X,Y, varname2, cmap='YlGnBu', transform=ccrs.PlateCarree(), vmin=0, vmax=600)
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
    plt.show()
-----------------------------------------------------------------------------------------------------------------------------
# Creates a 4 subplot figure when provided with the four variable names, lons, lats, and extent.
def four_fig(varname1, varname2, varname3, varname4, lons, lats, extent):
    fig, axes= plt.subplots(nrows=2, ncols=2, figsize=(25,20), subplot_kw={'projection':ccrs.PlateCarree()})
    extent = [-80, -50, 0, -35 ]
    level= np.linspace(-100,80,35)
    level1=np.linspace(0,800,40)

    ax1= axes[0][0]
    ax2= axes[0][1]
    ax3= axes[1][0]
    ax4= axes[1][1]

    ax1.add_feature(cfeature.LAND)
    ax1.add_feature(cfeature.COASTLINE, color='k')
    ax2.add_feature(cfeature.LAND)
    ax2.add_feature(cfeature.COASTLINE, color='k')
    ax1.set_extent(extent)
    ax2.set_extent(extent)
    ax3.set_title('MH-PI (JJA)', loc='left', fontsize=25)
    ax4.set_title('MH-PI (DJF)', loc='left', fontsize=25)
    ax3.add_feature(cfeature.LAND)
    ax3.add_feature(cfeature.COASTLINE, color='k')
    ax4.add_feature(cfeature.LAND)
    ax4.add_feature(cfeature.COASTLINE, color='k')
    ax3.set_extent(extent)
    ax4.set_extent(extent)
    ax1.set_title('PI (JJA)', loc='left', fontsize=25)
    ax2.set_title('PI (DJF)', loc='left', fontsize=25)

    X=lons
    Y=lats
    X, Y= np.meshgrid(X,Y)

    ax1.pcolor(X,Y, varname1, cmap= 'BrBG', transform=ccrs.PlateCarree())
    ax1.plot(-75.180, -14.5, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
    ax1.text(-75, -14, 'Palpa', fontsize= 20, transform=ccrs.PlateCarree())
    ax1.plot(-69.3354, -15.9254, marker= '.' , color= 'black', markersize= 15, transform= ccrs.PlateCarree())
    ax1.text(-69, -15, 'Titicaca', fontsize= 20, transform=ccrs.PlateCarree())
    ax2.pcolor(X,Y, varname2, cmap='YlGnBu', transform=ccrs.PlateCarree())
    ax3.pcolor(X,Y, varname3, cmap='BrBG', transform=ccrs.PlateCarree())
    ax4.pcolor(X,Y, varname4, cmap='BrBG', transform=ccrs.PlateCarree())
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

    #plt.title('Precipitation anomalies PI [JJA]')
    fig.canvas.draw()   # the only way to apply tight_layout to matplotlib and cartopy is to apply canvas firt 
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.88, top=0.94, bottom=0.10, hspace=0.1)
    plt.show()
------------------------------------------------------------------------------------------------------------------------------
# Creates a two_subplot image with quiver on top
def two_fig_quiver(varname1, varname2, u_vector, v_vector, u_vector_1, v_vector_1, lons, lats, extent):
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

    plot1=ax1.pcolormesh(X,Y, varname1, cmap='RdBu', transform=ccrs.PlateCarree())
    ax1.quiver(X[::2],Y[::2], u_vector[::2], v_vector[::2], transform= ccrs.PlateCarree() , color='white', scale=350)
    plot2=ax2.pcolormesh(X,Y, varname2, cmap='RdBu', transform=ccrs.PlateCarree())
    ax2.quiver(X[::2],Y[::2], u_vector_1[::2], v_vector_1[::2], transform= ccrs.PlateCarree() , color='white' , scale=100)
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




        
