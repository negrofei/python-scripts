"""
Mapa de la hora de máxima precipitación en todo el periodo.
Se calcula el ciclo medio diurno para cada punto de reticula cuando se supera algun umbral de pp
se extrae el tiempo de máxima precipitación. 

for execution:
    $ python3 map_time_of_maxPP_wholeperiod.py


"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Importo Paquetes

import os
import sys
import Docpy
from Docpy.functions import printer
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

# Agarro los rains de cada miembro
wrfs = [Docpy.WRF(path) for path in glob.glob('./rains*')] 

# Figs dir
fail = sys.argv[0][:-3]     # name of script
figsdir = os.path.abspath( os.path.join(wrfs[0].ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

# Otras cosas
acum = 3
lat, lon = wrfs[0].get_latlon()

# Spin up de 15 días 
spinup_time = int(15 * 24 / acum)
printer('Spin up time:', int(spinup_time*acum/24),'días')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Calculos
printer('Acum', acum)
precips = [Docpy.precip.calc_precip_acum(wrf, acum=acum)[spinup_time:-1] for wrf in wrfs]
# mm/day (saco el último  tiempo que son 0s)

### Mapa de la hora de la precipitación máxima
#    
#   * Necesito un mapa de colores circular
#   * Tengo que calcular el ciclo diurno promedio para cada punto de retícula
#   * Luego tengo que extraer la hora de la máxima precipitación de ese ciclo diurno promedio. 

# Ciclo diurno de precipitación mayor a un umbral en mm/3h
mean_cicles = [np.ma.zeros( (int(24/acum),)+precip.shape[1:] ) for precip in precips]

threshold = float(input('Threshold?\t'))

precip_threshold = [np.ma.masked_less(precip, threshold) for precip in precips]

for p, precip in enumerate(precips):
    for t in range(mean_cicles[p].shape[0]):
        mean_cicles[p][t] = np.ma.mean( precip_threshold[p][t::mean_cicles[p].shape[0]], axis=0)

# Promedio el ciclo diurno entre los miembros
mean_cicle = np.ma.mean(mean_cicles, axis=0)

# Maximo de cada punto de grilla
max_mean_cicle = np.ma.max(mean_cicle, axis=0)

# Agarro el tiempo en el que se cumple que llega al maximo
time_max_pp = np.ma.masked_equal(np.ones_like(max_mean_cicle),1) 

for t in range(0,mean_cicle.shape[0]):
    time_max_pp = np.ma.where( mean_cicle[t]==max_mean_cicle, acum*t, time_max_pp)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Plot Properties

plt.rcParams.update({'font.size': 14})
box = Docpy.functions.utils['boxHR3-44']    #Spatial spinup de 44 puntos (2 veces la relaxation)
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='k',
        facecolor='none',
        )
#provincias = Docpy.plot.mapping.load_provinces()
provincias = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_state_provinces_lines', 
        scale='50m', 
        facecolor='none', 
        edgecolor='gray')
# Niveles
levels = [0,3,6,9,12,15,18,21,24]
paleta = plt.get_cmap('twilight_shifted',len(levels)+1)
colores = [paleta(i) for i in range(paleta.N)]
colores[-1] = colores[-3]+np.diff(colores,axis=0)[-3]
new_paleta = colors.LinearSegmentedColormap.from_list('new_colormap', colores, paleta.N)
norm = colors.BoundaryNorm(levels, paleta.N)

### Plot
fig = plt.figure( figsize=(10,10) )
axes = plt.axes(projection=ccrs.PlateCarree())
mapa = axes.pcolormesh(
        lon, lat, time_max_pp,
        transform=ccrs.PlateCarree(),
        cmap=new_paleta,
        #extend='both',
        #levels=levels,
        norm=norm,
        )
axes.set_extent(box)
axes.add_feature(bordes)
#axes.add_feature(provincias)
axes.coastlines('50m', color='k')

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
cbar.set_ticks(np.array(levels)+1.5)
cbar.ax.yaxis.set_ticklabels([str(level) for level in levels]) 
cbar.set_label('UTC')

# Mas configs
axes.set_title(
        ' '.join(
            [
                'Time of max PP', 
                'above {} mm/3h'.format(threshold),
                os.path.split(wrfs[0].ruta)[-1],
                'acum{:02d}'.format(acum)
                ]
            )
        )

#plt.show()
sape = '_'.join(
        [
            'map_time_of_maxPP_OND',
            'acum{:02}'.format(acum),
            'above_{}mm3h'.format(threshold),
            wrfs[0].ruta.split('/')[-3],
            wrfs[0].ruta.split('/')[-2],
            ]
        )
fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
plt.close()

