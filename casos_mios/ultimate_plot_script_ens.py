#%% Importo Paquetes

import os
import sys
import numpy as np
import Docpy
import time
import re
##
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
##
from tqdm import tqdm
#%% Rutas y configs

fail = sys.argv[0].split('.')[0]        # nombre del script

#### rutas

ruta = os.path.abspath('../')
ruta_figs = os.path.join(ruta,'mean_of_configs', 'figuras')
try:
    os.mkdir(ruta_figs)
    print(ruta_figs, 'creado con exito')
except:
    print(ruta_figs, 'ya existe')

# definiciones
res_linestyle = {'2.4km':'dashed',
                 '4km':'solid',
                 '12km':'dashdot',
                 '20km':'dotted'}

simus = ['2.4km_anidado12',
         '4km_anidado12',
         '4km_anidado12_B',
         '4km_anidado20_B',
         '4km_B',
         '4km_B_15S47W',
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

size_colors = {'1A12':'tab:blue',
               '2B12':'tab:orange',
               '2B20':'chocolate',
               '2B'  :'goldenrod',
               '3':'tab:red'}
size_colors_names = dict(zip([*size_colors], ['1-n12','2-n12','2-n20','2-noNest','3-noNest']))

caso = ruta.split('/')[['caso2' in x for x in ruta.split('/')].index(True)]
datos = ['cmorph','mswep','imerg']
casito = caso[4:]

def load_wrfs(simus):
    size_colors_names = dict(zip([*size_colors], ['1-n12','2-n12','2-n20','2-noNest','3-noNest']))
    wrfs = {}
    for sim in simus:
        for dom in simus_d0s[sim]:
            ruta_wrf = os.path.abspath(
                    os.path.join(
                        ruta,
                        'mean_of_configs',
                        sim,
                        '_'.join(['wrfout','arw',dom,casito+'_boxed_ens.nc'])
                        )
                    )
            try:
                wrf = Docpy.WRF(ruta_wrf)
            except FileNotFoundError as err:
                print(err.strerror, ruta_wrf)
            else:
                namelist_path = '/home/martin.feijoo/casos_mios/'+caso+'/config_inicial/'+sim
                res = Docpy.functions.get_namelist_info(os.path.join(namelist_path, 'namelist.input'))[0]
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
                wrfs[sim,dom] = wrf
    return wrfs

wrfs = load_wrfs(simus)

while True:
    opt = input('BOX4 o BOX4B o HR3?? ([1]/2/3) \t')
    if (opt.lower() == '1') or (opt == ''):
        box = Docpy.functions.utils['box']
        box_name = 'box4_chica'
        break
    elif opt.lower() == '2':
        box = [-87.1313, -37.5094, -47.2198, -7.15]
        # esta caja es LR2
        box_name = 'box4_grande'
        break
    elif opt.lower() == '3':
        box = [-76.5693, -47.7057, -36.535, -15.3772]
        # esta caja es HR3
        box_name = 'box_HR3'
        break
    else:
        print('Escribí bien, pelotudo')

plt.rcParams.update({'font.size':14})
#############################%% Z850 (cont) + UV (vector) + V (shaded) ##############################
if False:
    for sim in simus:
        d = simus_d0s[sim][-1]
        print(sim)
        wrf = wrfs[sim,d]
        import scipy.ndimage as sci
        lat, lon = wrf.get_latlon()
        times, delta_t_hs = wrf.timedata() 
        every = 3
        cada = every//delta_t_hs
        # Vq
        Vq = Docpy.campos.calc_Vq_integrated(wrf, delta_t=every)
        # vientos 
        plevs = list(wrf.get_variables('lev'))
        #U = wrf.get_variables('umet')[::cada,plevs.index(850),:,:]
        V = wrf.get_variables('vmet')[::cada,plevs.index(850),:,:]
        vpaleta = plt.get_cmap('Spectral_r')
        vlevels = np.linspace(-30,30,9)
        # z850
        zlevels = np.arange(1400,1575,25)
        z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)[::cada,:,:]
        smooth = False
        if smooth == True:
            a = -8.75
            b = 1.55
            r = (np.diff(lon)[0]+np.diff(lat)[0])/2
            f = lambda x: a*x+b
            z850 = sci.filters.gaussian_filter(z850, f(r))
        z850 = np.ma.masked_greater(z850, 2000)
        # voy a poner las flechas donde las pone ERA-Interim
        if opt == '2':
            intervalo = 2
        else:
            intervalo = 1
        era = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/UV_ERAI_200503.nc')
        latera, lonera = era.get_latlon()
        imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box)
        lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
        lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])

        # voy a poner algunos contornos rojos para la precipitacion
        # aca tengo 2 opciones: o poner la precipitación acumulada cada 3 horas con un lag
        #                       o poner la precipitacion acumulada cada 6
        # voy a hacer las 2
        levels_pp = [20]
        for opt in ['acum_in_{cada}'.format(cada=every),'acum_in_running']:
            if opt == 'acum_in_{cada}'.format(cada=every):
                precip = Docpy.precip.calc_precip_acum(wrf, acum=every)
                lag = 1
            elif opt == 'acum_in_running':
                precip = Docpy.precip.calc_running_acum(wrf, acum=every, run=every*2)
                lag = 0
            for t in tqdm(range(len(times[::cada])-lag-1), desc='V + z850'):
                fig = plt.figure(figsize=(10,10))
                axes = plt.axes(projection=ccrs.PlateCarree())
                axes.set_extent(box)
                mapa = axes.contourf(lon, lat, V[t],
                                    transform=ccrs.PlateCarree(),
                                    cmap=vpaleta,
                                    extend='both',
                                    levels=vlevels)
                n = 200
                flechas = axes.quiver(lon[lista_lon[::intervalo]], lat[lista_lat[::intervalo]], 
                        Vq[0][:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                        Vq[1][:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                                      transform=ccrs.PlateCarree(),
                                      scale=5, scale_units='inches')
                axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{kg}{m s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')

                contorno = axes.contour(lon, lat, z850[t],
                                        transform=ccrs.PlateCarree(),
                                        colors='k',
                                        levels=zlevels,
                                        )
                axes.clabel(contorno, inline=True, fmt='%0d')

                contorno2 = axes.contour(lon, lat, precip[t+lag,:,:],
                                        transform=ccrs.PlateCarree(),
                                        colors='purple',
                                        levels=levels_pp,
                                        alpha=0.8,
                                        )
                axes.clabel(contorno2, inline=True, fmt='%0d')
                bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                        category='cultural',
                        name='admin_0_boundary_lines_land',
                        scale='50m',
                        edgecolor='gray',
                        facecolor='none'
                        )
                axes.add_feature( bordes )
                axes.coastlines('50m', color='gray')
                #for rec in provincias.records():
                #    axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
                
                #plt.show()
                sape = '_'.join(['Vq_cfhi',opt, box_name, wrf.dominio, caso, sim, times[::int(every/delta_t_hs)][t]])
                fig.savefig(os.path.join(ruta_figs,sape+'.png'), dpi=300, bbox_inches='tight')
                plt.close()

#### z850 y vientos en vectores ####
if True:
    for sim in simus[1:]:
        d = 'd01'
        print(sim)
        wrf = wrfs[sim,d]
        import scipy.ndimage as sci
        lat, lon = wrf.get_latlon()
        times, delta_t_hs = wrf.timedata()
        every = 3
        cada = every//delta_t_hs

        # vientos 
        plevs = list(wrf.get_variables('lev'))
        U = wrf.get_variables('umet')[::cada,plevs.index(850),:,:]
        V = wrf.get_variables('vmet')[::cada,plevs.index(850),:,:]
        # z850
        zlevels = np.arange(1400,1575,25)
        z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)[::cada,:,:]
        z850 = np.ma.masked_greater(z850, 2000)
        zpaleta = plt.get_cmap('cool')
        # voy a poner las flechas donde las pone ERA-Interim
        if opt == '2':
            intervalo = 2
        else:
            intervalo = 1
        era = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/UV_ERAI_200503.nc')
        latera, lonera = era.get_latlon()
        imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box)
        lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
        lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])

        for t in tqdm(range(len(times[::cada])), desc='V + z850'):
            fig = plt.figure(figsize=(10,10))
            axes = plt.axes(projection=ccrs.PlateCarree())
            axes.set_extent(box)
            mapa = axes.contourf(lon, lat, z850[t],
                                transform=ccrs.PlateCarree(),
                                cmap=zpaleta,
                                extend='both',
                                levels=zlevels)
            n = 15
            flechas = axes.quiver(lon[lista_lon[::intervalo]], lat[lista_lat[::intervalo]],
                    U[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                    V[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                                  transform=ccrs.PlateCarree(),
                                  scale=5, scale_units='inches')
            axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')

            contorno = axes.contour(lon, lat, z850[t],
                                    transform=ccrs.PlateCarree(),
                                    colors='yellow',
                                    levels=zlevels,
                                    alpha=0.8,
                                    )
            #axes.clabel(contorno, inline=True, fmt='%0d')

            bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                    category='cultural',
                    name='admin_0_boundary_lines_land',
                    scale='50m',
                    edgecolor='darkgray',
                    facecolor='none'
                    )
            axes.add_feature( bordes )
            axes.coastlines('50m', color='gray')
            #for rec in provincias.records():
            #    axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
            cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=zlevels)#, format=OOMFormatter(-3))
            cbar.set_label(r'm')

            if t == 0:
                plt.show()
            sape = '_'.join(['z850_UV',opt, box_name, wrf.dominio, caso, sim, times[::cada][t]])
            fig.savefig(os.path.join('.',sape+'.png'), dpi=300, bbox_inches='tight')
            plt.close()

