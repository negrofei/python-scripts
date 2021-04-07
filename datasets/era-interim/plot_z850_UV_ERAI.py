"""
script que grafica ERAI para el trabajo de dominios
como el script 'ultimate_plot_script.py' 
"""

import os
import sys
import numpy as np
import Docpy
import re
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from tqdm import tqdm 
import scipy.ndimage as sci

fail = sys.argv[0].split('.py')[0]
ruta_figs = os.path.abspath( os.path.join('..','figuras',fail) )
if not os.path.exists(ruta_figs):
    os.makedirs(ruta_figs)
# cargo ERAI
era = Docpy.ERAI(sys.argv[1])
#era_winds = Docpy.ERAI(sys.argv[2])
lat, lon = era.get_latlon()

#props
caso = re.search(r'\d{6}',era.ruta).group()
plt.rcParams.update({'font.size':14})

# graph
while True:
    opt = input('BOX4  o BOX4B?? ([1]/2) \t ')
    if (opt.lower() == '1') or (opt==''):
        box = Docpy.functions.utils['box']
        box_name = 'box4_chica'
        intervalo = 1
        break
    elif opt.lower() == '2':
        box = [-87.1313, -37.5094, -47.2198, -7.15]
        # esta caja es LR2
        box_name = 'box4_grande'
        intervalo = 2
        break
    else:
        print('Escribí bien, pelotudo')
times, delta_t_hs = era.timedata()
every = 6
cada = int(every//delta_t_hs)
# Vq
Vq = Docpy.campos.calc_Vq_integrated(era, delta_t=every)

# vientos
plevs = list(era.get_levels())
V = era.get_variables(era.dvars_era['V'])[::cada,plevs.index(85000),:,:]
vpaleta = plt.get_cmap('Spectral_r')
vlevels = np.linspace(-30,30,9)

# z850
zlevels = np.arange(1400,1575,25)
lonmin_idx = np.where( np.abs(lon-box[0]) == np.min(np.abs(lon-box[0])) )[0][0]
lonmax_idx = np.where( np.abs(lon-box[1]) == np.min(np.abs(lon-box[1])) )[0][0]
latmin_idx = np.where( np.abs(lat-box[2]) == np.min(np.abs(lat-box[2])) )[0][0]
latmax_idx = np.where( np.abs(lat-box[3]) == np.min(np.abs(lat-box[3])) )[0][0]
z850 = Docpy.campos.calc_z_lvl(era, lvl=85000)[::cada,latmin_idx:latmax_idx+1,lonmin_idx:lonmax_idx+1]
smooth = True
if smooth == True:
    a = -8.75
    b = 1.55
    r = (np.diff(lon)[0]+np.diff(lat)[0])/2
    f = lambda x: a*x+b
    z850 = sci.filters.gaussian_filter(z850, f(r))

# voy a poner las flechas donde las pone ERA-Interim
erai = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/UV_ERAI_200503.nc')
latera, lonera = erai.get_latlon()
imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box)
lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])

# plot Vq, V y Z
for t in tqdm(range(len(times[::cada])), desc='V + z850'):
    fig = plt.figure(figsize=(10,10))
    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.set_extent(box)
    mapa = axes.contourf(lon, lat, V[t],
            transform=ccrs.PlateCarree(),
            cmap=vpaleta,
            extend='both',
            levels=vlevels,
            )
    n = 200
    flechas = axes.quiver(lon[lista_lon[::intervalo]], lat[lista_lat[::intervalo]],
            Vq[0][:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
            Vq[1][:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                          transform=ccrs.PlateCarree(),
                          scale=5, scale_units='inches',
                          )
    axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{kg}{m s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')
    contorno = axes.contour(lon[lonmin_idx:lonmax_idx+1], lat[latmin_idx:latmax_idx+1], z850[t],
                            transform=ccrs.PlateCarree(),
                            colors='k',
                            levels=zlevels,
                            )
    axes.clabel(contorno, inline=True, fmt='%0d')
    bordes = cfeature.NaturalEarthFeature(          # limites de los paises
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale='50m',
            edgecolor='gray',
            facecolor='none',
            )
    axes.add_feature( bordes )
    axes.coastlines('50m', color='gray')
    
    gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 14, 'rotation': 25}
    gl.ylabel_style = {'size': 14, 'rotation': 25}
    axes.set_title(times[::int(every/delta_t_hs)][t], fontsize=14)

    # meto el cbar a mano
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=vlevels)#, format=OOMFormatter(-3))
    cbar.set_label(r'm s$^{-1}$')
    
    if t == 10:
        plt.show()
    sape = '_'.join(['Vq_cfhi', era.name, caso, box_name, times[::int(every/delta_t_hs)][t]])
    fig.savefig(os.path.join(ruta_figs,sape+'.png'), dpi=300, bbox_inches='tight')
    plt.close()


######### z850 con vientos #########
U = era.get_variables(era.dvars_era['U'])[::cada,plevs.index(85000),:,:]
V = era.get_variables(era.dvars_era['V'])[::cada,plevs.index(85000),:,:]
paleta = plt.get_cmap('cool')
for t in tqdm(range(0,len(times),cada), desc='Z850'):
    fig = plt.figure(figsize=(10,10))
    axes = plt.axes(projection=ccrs.PlateCarree())
    mapa = axes.contourf(lon[lonmin_idx:lonmax_idx+1], lat[latmin_idx:latmax_idx+1], z850[t,:,:],
            transform=ccrs.PlateCarree(),
            cmap=paleta,
            extend='both',
            levels=zlevels)

    n = 15
    flechas = axes.quiver(lon[lista_lon[::intervalo]], lat[lista_lat[::intervalo]],
            U[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
            V[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
            transform=ccrs.PlateCarree(),
            scale=5, scale_units='inches')
    axes.quiverkey(flechas,.875,1.025,1,'{}'.format(n)+r'$\frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')

    contorno = axes.contour(lon[lonmin_idx:lonmax_idx+1], lat[latmin_idx:latmax_idx+1], z850[t,:,:],
            transform=ccrs.PlateCarree(),
            levels=zlevels, color='k')
    axes.set_extent(box)
    bordes = cfeature.NaturalEarthFeature(          # limites de los paises
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale='50m',
            edgecolor='gray',
            facecolor='none'
            )
    axes.add_feature( bordes )
    axes.coastlines('50m', color='gray')
    #        for rec in provincias.records():
    #            axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
    gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 14, 'rotation': 25}
    gl.ylabel_style = {'size': 14, 'rotation': 25}
    
    axes.set_title(times[t], fontsize=14)
    
    # meto el cbar a mano
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=zlevels)
    cbar.set_label('m')
    
    #plt.show()
    sape = '_'.join(['z850', era.name, caso, box_name,times[t]])
    fig.savefig(os.path.join(ruta_figs,sape+'.png'), dpi=300, bbox_inches='tight')
    plt.close()

                                                                            
