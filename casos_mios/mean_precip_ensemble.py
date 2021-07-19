"""
Para cada simulacion promedio la precipitation de todas las configuraciones y grafico el mapa de precipitation.  
"""


#%% Importo Paquetes

import os
import sys
import numpy as np
import Docpy
from Docpy import WRF, Dato
from Docpy.plot.mapping import map_cbar
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

fail = sys.argv[0].split('.')[0]            #nombre del script

ruta = os.path.abspath('../')
ruta_figs = os.path.join(ruta, 'figuras')
try:
    os.mkdir(ruta_figs)
    print(ruta_figs, 'creado con exito')
except:
    print(ruta_figs, 'ya existe')

# props
caso = re.search(r'(caso\d{4}[-]\d{2})',ruta).group()
casito = caso[4:]
plt.rcParams.update({'font.size': 20})
acum_list = [3]#, 6, 24]                      # grafico para acumular cada 3 6 y 24

simus = [
        '2.4km_anidado12',
        '4km_anidado12',
        '4km_anidado12_B',
        '4km_anidado20_B',
        '4km_B',
        '4km_B_15S47W',
        ]
configs = [
        'config_inicial',
        'config_inicialYSU',
        'igual1_con-shaYSU',
        'pbl_MYNN',
        'shallow-on_MYNN_grell',
        'shallow-on_MYNN_grell_MP',
        ]

simus_d0s = dict(zip(simus, [['d02'],['d01','d02'],['d01','d02'],['d01','d02'],['d01'],['d01']]))

domains = dict(zip(simus,['1A12','1A12','2B12','2B20','2B','3']))

size_colors = {'1A12':'tab:blue',
        '2B12':'tab:orange',
        '2B20':'chocolate',
        '2B'  :'goldenrod',
        '3':'tab:red'}
size_colors_names = dict(zip([*size_colors], ['1(12)','2(12)','2(20)','2','3']))


wrfs = {}                                   # cargo todas las simulaciones
for sim in simus:
    Docpy.functions.printer(sim)
    for dom in simus_d0s[sim]:
        for config in configs:
            ruta_wrf = os.path.abspath(
                    os.path.join(
                        ruta,
                        '..',
                        config,
                        sim,
                        '_'.join(['wrfout','arw',dom,casito+'_boxed.nc'])
                        )
                    )
            try:
                wrf = WRF(ruta_wrf)
            except FileNotFoundError as err:
                print(err.strerror,ruta_wrf)
            else:
                wrf.config = config
                orig_path = os.path.abspath(wrf.ruta.split('interpolated_files')[0]+wrf.ruta.split('interpolated_files')[1])
                res = Docpy.functions.get_namelist_info(os.path.join(orig_path,'namelist.input'))[0]
                wrf.color = size_colors[domains[sim]]
                r = res[int(dom[-1])-1]
                if r  > 4000:
                    wrf.size = 'DR'
                else:
                    wrf.size = 'CP'
                wrf.sim = sim
                if r == 2400:
                    wrf.label = wrf.size+size_colors_names[domains[sim]][0]+'-'+str(r/1000)+size_colors_names[domains[sim]][1:]
                else:
                    wrf.label = wrf.size+size_colors_names[domains[sim]][0]+'-'+str(int(r/1000))+size_colors_names[domains[sim]][1:]
                if wrf.size == 'DR':
                    wrf.label = wrf.label[:-4]
                wrfs[config,sim,dom] = wrf

box = Docpy.functions.common_doms([os.path.join(wrfs[wrf].ruta,wrfs[wrf].filename) for wrf in [*wrfs]])

paleta = plt.get_cmap('gist_ncar')
paleta.set_under('white')

# voy a armar el ensemble para cada dominio
# primero tengo que hacer esto para cada acumulacion
for acum in acum_list:
    levels = Docpy.precip.d_format['levels'][acum]
    for sim in simus:
        Docpy.functions.printer(sim)
        ruta_figs_sims = os.path.join(ruta_figs, sim)
        try:
            os.makedirs( ruta_figs_sims )
        except FileExistsError as err:
            print(err.strerror,ruta_figs_sims, 'already exists')
        for dom in simus_d0s[sim]:
            lat, lon = wrfs['config_inicial',sim, dom].get_latlon()
            times = Docpy.functions.acum_time(wrfs['config_inicial', sim, dom], acum)
            precip_sim = np.zeros([len(configs), len(times), len(lat), len(lon)])
            for c,config in enumerate(configs):
                wrf = wrfs[config, sim, dom]
                precip_sim[c,:,:,:] = Docpy.precip.calc_precip_acum(wrf, acum)
            mean_precip = np.mean(precip_sim, axis=0)
            for t in tqdm(range(len(times)), desc='{:02}'.format(acum)):
                fig = plt.figure(figsize=(10,10))
                axes = plt.axes(projection=ccrs.PlateCarree())
                try:
                    mapa = axes.contourf(lon, lat, mean_precip[t,:,:],
                            transform=ccrs.PlateCarree(),
                            cmap=paleta,
                            extend='both',
                            levels=levels,
                            )
                except ValueError as err:
                    print('Hay un bardo con tiempo',times[t],'de la simulacion',sim+dom)
                axes.set_extent(box)
                bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                        category='cultural',
                        name='admin_0_boundary_lines_land',
                        scale='50m',
                        edgecolor='k',
                        facecolor='none',
                        )
                axes.add_feature( bordes )
                axes.coastlines('50m')
                gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style = {'size': 20, 'rotation': 25}
                gl.ylabel_style = {'size': 20, 'rotation': 25}
                
                axes.set_title(wrf.label, fontsize=30)
                
                # meto el cbar a mano
                #fig.subplots_adjust(right=0.8)
                #cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
                #cbar = fig.colorbar(mapa, cax=map_cbar(fig, axes, width=0.03), ticks=levels)
                #cbar.set_label('mm')

                # guardo la figura en un directorio adentro del directorio ruta_figs con el nombre de la simulacion
                fmts = ['.png','.pdf']
                for fmt in fmts:
                    sape = '_'.join(['ensemble','acum{:02}'.format(acum), wrf.dominio, caso, sim, times[t]])
                    fig.savefig(os.path.join(ruta_figs_sims,sape+fmt), dpi=300, bbox_inches='tight')
                plt.close()


            ### GUARDO EL CBAR APARTE ###
            fig, ax = plt.subplots(figsize=(10,10))
            mapa = ax.contourf(lon, lat, mean_precip[t,...], cmap=paleta, extend='both', levels=levels)
            cbar = fig.colorbar(mapa, cax=map_cbar(fig, ax, width=0.03), ticks=levels)
            cbar.set_label('mm')
            ax.remove()
            for fmt in fmts:
                sape = '_'.join(['ensemble_cbar','acum{:02}'.format(acum), wrf.dominio, caso, sim])
                fig.savefig(os.path.join(ruta_figs_sims,sape+fmt), dpi=300, bbox_inches='tight')

