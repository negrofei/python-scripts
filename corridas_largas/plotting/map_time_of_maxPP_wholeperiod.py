"""
Mapa de la hora de máxima precipitación en todo el periodo.
Se calcula el ciclo medio diurno para cada punto de reticula
se extrae el tiempo de máxima precipitación. 

for execution:
    $ python3 map_time_of_maxPP_wholeperiod.py


"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Importo Paquetes

import os
import sys
import Docpy
import time
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.colors as colors
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from tqdm import tqdm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
### Rutas y configs

wrfs = [Docpy.WRF(path) for path in glob.glob('./rains*')] 
# rains_acum03_alltimes_wrfout_<dir_of_chunks>.nc
acum = 3
lat, lon = wrfs[0].get_latlon()
fail = sys.argv[0][:-3]     # name of script
figsdir = os.path.abspath( os.path.join(wrfs[0].ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Calculos
Docpy.functions.printer(acum)
precips = [Docpy.precip.calc_precip_acum(wrf, acum=acum)[:-1] for wrf in wrfs]
# mm/day (saco el último  tiempo que son 0s)

### Mapa de la hora de la precipitación máxima
#    
#   * Necesito un mapa de colores circular
#   * Tengo que calcular el ciclo diurno promedio para cada punto de retícula
#   * Luego tengo que extraer la hora de la máxima precipitación de ese ciclo diurno promedio. 


# Ciclo diurno
mean_cicles = [np.zeros( (int(24/acum),)+precip.shape[1:] ) for precip in precips]

for p, precip in enumerate(precips):
    for t in range(mean_cicles[p].shape[0]):
        mean_cicles[p][t] = np.mean( precip[t::mean_cicles[p].shape[0]], axis=0)

# Calculo el maximo del ciclo diurno
max_mean_cicles = [np.max(mean_cicle, axis=0) for mean_cicle in mean_cicles]

# Agarro el tiempo en el que se cumple que llega al maximo
time_max_pps = [np.zeros_like(max_mean_cicle) for max_mean_cicle in max_mean_cicles]

for m, mean_cicle in enumerate(mean_cicles):
    for t in range(1,mean_cicle.shape[0]):
        time_max_pps[m] = np.where( mean_cicle[t]==max_mean_cicles[m], acum*t, time_max_pps[m])

# Promedio en caso que sea un ensemble

time_max_pp = np.mean(
        np.stack(
            time_max_pps,
            axis=0,
            ),
        axis=0,
        )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Plot Properties

plt.rcParams.update({'font.size': 14})
box = Docpy.functions.utils['boxHR3']
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='k',
        facecolor='none',
        )
paleta = plt.get_cmap('twilight_shifted')

# Niveles
levels = [0,3,6,9,12,15,18,21]

### Plot
fig = plt.figure( figsize=(10,10) )
axes = plt.axes(projection=ccrs.PlateCarree())
mapa = axes.contourf(
        lon, lat, time_max_pp,
        transform=ccrs.PlateCarree(),
        cmap=paleta,
        extend='both',
        levels=levels,
        )
axes.set_extent(box)
axes.add_feature(bordes)
axes.coastlines('50m', color='gray')

# Gridlines
gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14, 'rotation': 25}
gl.ylabel_style = {'size': 14, 'rotation': 25}

# Meto el cbar a mano
fig.subplots_adjust(right=0.8)
cbar_xmin = axes.get_position().xmax+0.05
cbar_ymin = axes.get_position().ymin
cbar_width = 0.025
cbar_hight = axes.get_position().ymax-axes.get_position().ymin
cbar_ax = fig.add_axes([cbar_xmin, cbar_ymin, cbar_width, cbar_hight])
cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
cbar.set_label('UTC')

# Mas configs
axes.set_title(' '.join(['Time of max PP', os.path.split(wrfs[0].ruta)[-1],'acum{:02d}'.format(acum)]))

#plt.show()
sape = '_'.join(
        [
            'map_time_of_maxPP_OND',
            'acum{:02}'.format(acum),
            wrfs[0].ruta.split('/')[-3],
            wrfs[0].ruta.split('/')[-2],
            ]
        )
fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
plt.close()

