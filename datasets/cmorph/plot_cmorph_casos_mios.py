#!/usr/bin/env python3

"""
Script to plot CMORPH precipitation given a 'acum' acumulation parameter

for execution:
    $ python3 plot_cmorph.py <CMORPH_dataset>


"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Importo Paquetes

import os
import sys
import numpy as np
import Docpy 
import time
from tqdm import tqdm

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Inputs

filename = sys.argv[1]      # name of the file
fail = sys.argv[0][:-3]     # name of script
ruta_figs = os.path.join( os.path.split(os.path.abspath(filename))[0],'figuras', fail)
# me quedo con la raiz de la ruta del archivo 
try:
    os.makedirs(ruta_figs)
except:
    print(ruta_figs,'ya existe')

acum_list = [3, 6, 24]

### Read file
cmorph = Docpy.Dato(os.path.abspath(filename))


### Plot Properties and utils
lat, lon = cmorph.get_latlon()
plt.rcParams.update({'font.size': 14})
paleta = plt.get_cmap('gist_ncar')
paleta.set_under('white')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Plot
Docpy.functions.printer('Grafico precipitation...')

for acum in acum_list:
    levels = Docpy.precip.d_format['levels'][acum]
    times = Docpy.functions.acum_time(cmorph, acum=acum)
    if '200503' in filename:
        tini = times.index('2005-03-11 00:00:00')
        tfin = times.index('2005-03-14 00:00:00')
    if '201511' in filename:
        tini = times.index('2015-11-09 00:00:00')
        tfin = times.index('2015-11-12 00:00:00')
    if '201610' in filename:
        tini = times.index('2016-10-23 00:00:00')
        tfin = times.index('2016-10-26 00:00:00')
    times =  times[tini:tfin]
    precip = Docpy.precip.calc_precip_acum(cmorph, acum)[tini:tfin,:]
    for t in tqdm(range(len(times)), desc='{:02}'.format(acum)):
        fig = plt.figure(figsize=(10,10))
        axes = plt.axes(projection=ccrs.PlateCarree())
        mapa = axes.contourf(lon, lat, precip[t,:,:],
                             transform=ccrs.PlateCarree(),
                             cmap=paleta,
                             extend='both',
                             levels=levels,
                             )

        axes.set_extent(Docpy.functions.utils['box'])
        # paises
        bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='50m',
                edgecolor='k',
                facecolor='none'
                )
        axes.add_feature( bordes )
        axes.coastlines(resolution='50m')
#        for rec in provincias.records():
#            axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
        gl = axes.gridlines(draw_labels=True,color='black',alpha=0.2,linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14, 'rotation':25}
        gl.ylabel_style = {'size': 14}
        
        axes.set_title(times[t], fontsize=14)

        # meto el cbar a mano
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
        cbar.set_label('mm')
        
        #plt.show()
        sape = '_'.join([this_file,'cmorph','acum{:02d}'.format(acum),times[t][:-6]])
        fig.savefig(os.path.join(save_dir,sape+'.jpg'), dpi=300, bbox_inches='tight')
        plt.close()

            
Docpy.functions.printer('Terminado')
