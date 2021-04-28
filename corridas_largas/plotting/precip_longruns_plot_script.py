"""
Grafico mapas de todos los tiempos de la precipitation de la corrida larga

for execution:

    $ python3 precip_longruns_plot_script.py rains_acum03_alltimes_wrfout_<dir_of_chunks>.nc



"""


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Importo Paquetes

import os
import sys
import numpy as np
import Docpy
import time

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from tqdm import tqdm


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
### Rutas y configs

wrf = Docpy.WRF(sys.argv[1])    # rains_acum03_alltimes_wrfout_<dir_of_chunks>.nc
acum_list = [3, 6, 24]          # grafico para acumular cada 3 6 y 24
lat, lon = wrf.get_latlon()
fail = sys.argv[0][:-3]         # name of scripts 
figsdir = os.path.abspath( os.path.join(wrf.ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

### Plot Properties

plt.rcParams.update({'font.size': 14})
box = Docpy.functions.utils['boxHR3']
paleta = plt.get_cmap('gist_ncar')
paleta.set_under('white')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Plot
for acum in acum_list:
    levels = Docpy.precip.d_format['levels'][acum]
    precip = Docpy.precip.calc_precip_acum(wrf, acum=acum)
    times = Docpy.functions.acum_time(wrf, acum=acum)
    for t in tqdm(range(len(times)), desc='{:02}'.format(acum)):
        fig = plt.figure(figsize=(10,10))
        axes = plt.axes(projection=ccrs.PlateCarree())
        mapa = axes.contourf(
                lon, lat, precip[t,:,:],
                transform=ccrs.PlateCarree(),
                cmap=paleta,
                extend='both',
                levels=levels,
                )
        # Limits of plot
        axes.set_extent(box)
        
        # Countries borders
        bordes = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='50m',
                edgecolor='gray',
                facecolor='none',
                )
        axes.add_feature( bordes )
        axes.coastlines('50m', color='gray')

        # Gridlines
        gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14, 'rotation': 25}
        gl.ylabel_style = {'size': 14, 'rotation': 25}

        # Title
        axes.set_title(times[t], fontsize=14)

        # meto el cbar a mano
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
        cbar.set_label('mm')

        #plt.show()
        sape = '_'.join(
                [
                    'precip',
                    'acum{:02}'.format(acum),
                    wrf.dominio,
                    wrf.ruta.split('/')[-3],
                    wrf.ruta.split('/')[-2],
                    times[t],
                    ]
                )
        fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
        plt.close()

