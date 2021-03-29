import Docpy
import matplotlib.pyplot as plt
import sys
import os

import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
# current file

fail = sys.argv[0].split('.py')[0]

# rutas

ruta = os.path.abspath('.')
ruta_figs = os.path.join(ruta,'../figuras')
try:
    os.mkdir(ruta_figs)
    print(ruta_figs, 'creado con éxito')
except:
    print(ruta_figs, 'ya existe')

# definiciones
res_linestyle = {
        '2.4km':'dashed',
        '4km':'solid',
        '12km':'dashdot',
        '20km':'dotted'
        }

simus = [
        #'2.4km_anidado12',
        #'4km_anidado12',
        #'4km_anidado12_B',
        '4km_anidado20_B',
        '4km_B',
        '4km_B_15S47W'
        ]

configs = [
        'config_inicial',
        #'config_inicialYSU',
        #'igual1_con-shaYSU',
        'pbl_MYNN',
        #'shallow-on_MYNN_grell',
        'shallow-on_MYNN_grell_MP',
        ]


simus_d0s = {
        '2.4km_anidado12':['d02'],
        '4km_anidado12':['d01','d02'],
        '4km_anidado12_B':['d01','d02'],
        '4km_anidado20_B':['d01','d02'],
        '4km_B':['d01'],
        '4km_B_15S47W':['d01'],
        }

domains = dict(zip([*simus_d0s],['1A12','1A12','2B12','2B20','2B','3']))

caso = ruta.split('/')[['caso2' in x for x in ruta.split('/')].index(True)]
datos = ['era5','erai']
casito = caso[4:]

# Funcion para cargar todas las simulaciones para cada config
def load_wrfs(simus, configs):
    size_colors = {
            '1A12':'tab:blue',
            '2B12':'tab:orange',
            '2B20':'chocolate',
            '2B'  :'goldenrod',
            '3':'tab:red'
            }
    size_colors_names = dict(zip([*size_colors], ['1-n12','2-n12','2-n20','2-noNest','3-noNest']))
    wrfs = {}
    for sim in simus:
        for dom in simus_d0s[sim]:
            for config in configs:
                ruta_wrf = os.path.abspath(
                        os.path.join(
                            ruta,
                            '../../',
                            config,
                            sim,
                            '_'.join(['wrfout','arw',dom,casito+'.nc'])
                            )
                        )
                try:
                    wrf = Docpy.WRF(ruta_wrf)
                except FileNotFoundError as err:
                    print(err.strerror, ruta_wrf)
                else:
                    wrf.config = config
                    res = Docpy.functions.get_namelist_info(os.path.join(wrf.ruta, 'namelist.input'))[0]
                    wrf.color = size_colors[domains[sim]]
                    r = res[int(dom[-1])-1]
                    if r == 2400:
                        wrf.linestyle = res_linestyle[str(r/1000)+'km']
                    else:
                        wrf.linestyle = res_linestyle[str(int(r/1000))+'km']
                        
                    if r  > 4000:
                        wrf.size = 'LR'
                    else:
                        wrf.size = 'HR'
                    wrf.sim = sim
                    if r == 2400:
                        wrf.label = str(r/1000)+'km'+'_'+wrf.size+size_colors_names[domains[sim]]
                    else:
                        wrf.label = str(int(r/1000))+'km'+'_'+wrf.size+size_colors_names[domains[sim]]
                    wrfs[config,sim,dom] = wrf
    return wrfs

wrfs = load_wrfs(simus, configs)

# niveles verticales
lvls = list(wrfs[[*wrfs][0]].get_levels())

# Defino en qué region voy a calcular el promedio de la temperatura

box = Docpy.functions.common_doms([os.path.join(wrfs[wrf].ruta,wrfs[wrf].filename) for wrf in [*wrfs]]) # pruebo primero calcularlo en toda la region comun

target_box = Docpy.functions.utils['box']

# calculo la lista de tiempos
time_wrf = []
for wrf in [*wrfs]:
    time_wrf.append( wrfs[wrf].acum_time(3) )
time_common = Docpy.functions.common_times(time_wrf)
plt.rcParams.update({'font.size':14})

# calculo el espesor 1000 - 500 
for k in [*wrfs]:
    print(k)
    wrf = wrfs[k]
    lat, lon = wrf.get_latlon(box)
    tiempos, delta_t_hs = wrf.timedata()
    inter = int(3/delta_t_hs)   #cada 3 hs hago esto
    z = wrf.get_variables('geopt', box)[::inter,...]/9.81
    z850  = z[:,lvls.index(850),...]
    delta_z = z[:,lvls.index(500),...] - z[:,lvls.index(950),...]
    for t in range(z.shape[0]):
        fig = plt.figure(figsize=(10,7))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(box)
        shade = ax.contourf(lon, lat, delta_z[t], transform=ccrs.PlateCarree(), levels=np.linspace(5130,5410,8), extend='both')
        contour = plt.contour(lon, lat, z850[t],transform=ccrs.PlateCarree(), levels=np.arange(1400,1575,25), colors='k')

        bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='50m',
                edgecolor='k',
                facecolor='none'
                )
        ax.add_feature(bordes)
        ax.coastlines('50m')
        # dibujo el target domain
        ax.hlines(y=target_box[2], xmin=target_box[0], xmax=target_box[1], transform=ccrs.PlateCarree())
        ax.hlines(y=target_box[3], xmin=target_box[0], xmax=target_box[1], transform=ccrs.PlateCarree())
        ax.vlines(x=target_box[0], ymin=target_box[2], ymax=target_box[3], transform=ccrs.PlateCarree())
        ax.vlines(x=target_box[1], ymin=target_box[2], ymax=target_box[3], transform=ccrs.PlateCarree())

        ax.clabel(contour, inline=True, fmt='%0d')
        ax.set_title(wrf.label+' '+tiempos[::inter][t])
        
        # meto el cbar a mano
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85,0.15, 0.02, 0.7])
        cbar = fig.colorbar(shade, cax=cbar_ax, ticks=shade.levels)
        cbar.set_label('m')

        fig.suptitle(wrf.config)
        fig.savefig(os.path.join(ruta_figs,
            '_'.join(['espesor',wrf.label,wrf.config,fail, caso, tiempos[::inter][t]]))+'.png', dpi=150, bbox_inches='tight')
        if t%5==0:
            pass
            #plt.show()
        plt.close()

# Cargo los datos de ERA5 y ERAI

ruta_erai = '/home/martin.feijoo/ERA-Interim/casos_mios/'

# ERA5
ruta_era5 = '/home/martin.feijoo/ERA5/casos_mios/ERA5_'+casito+'_UVqZ.nc'
era5 = Docpy.ERA5(ruta_era5)
latera, lonera = era5.get_latlon(box)
lvlsera5 = list(era5.get_levels())
tiempo_era5, delta_t_hs_era5 = era5.timedata()
inter = int(3/delta_t_hs_era5)   #cada 3 hs hago esto
z_era5 = era5.get_variables('z', box)[::inter,...]/9.81
z850_era5 = z_era5[:,lvlsera5.index(850),...]
delta_z_era5 = z_era5[:,lvlsera5.index(500),...] - z_era5[:,lvlsera5.index(950),...]

for t in range(z_era5.shape[0]):
    fig = plt.figure(figsize=(10,7))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(box)
    shade = ax.contourf(lonera, latera, delta_z_era5[t], transform=ccrs.PlateCarree(), levels=np.linspace(5130,5410,8), extend='both')
    contour = plt.contour(lonera, latera, z850_era5[t],transform=ccrs.PlateCarree(), levels=np.arange(1400,1575,25), colors='k')

    bordes = cfeature.NaturalEarthFeature(          # limites de los paises
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale='50m',
            edgecolor='k',
            facecolor='none'
            )
    ax.add_feature(bordes)
    ax.coastlines('50m')
    # dibujo el target domain
    ax.hlines(y=target_box[2], xmin=target_box[0], xmax=target_box[1], transform=ccrs.PlateCarree())
    ax.hlines(y=target_box[3], xmin=target_box[0], xmax=target_box[1], transform=ccrs.PlateCarree())
    ax.vlines(x=target_box[0], ymin=target_box[2], ymax=target_box[3], transform=ccrs.PlateCarree())
    ax.vlines(x=target_box[1], ymin=target_box[2], ymax=target_box[3], transform=ccrs.PlateCarree())

    ax.clabel(contour, inline=True, fmt='%0d')
    ax.set_title('ERA5'+' '+tiempo_era5[::inter][t])
    
      # meto el cbar a mano
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85,0.15, 0.02, 0.7])
    cbar = fig.colorbar(shade, cax=cbar_ax, ticks=shade.levels)
    cbar.set_label('m')

    fig.savefig(os.path.join(ruta_figs,
        '_'.join(['espesor','ERA5',fail, caso, tiempo_era5[::inter][t]]))+'.png', dpi=150, bbox_inches='tight')
    if t%5==0:
        pass
        #plt.show()
    plt.close()
