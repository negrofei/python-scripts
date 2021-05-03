"""
Mapa de la precipitación media de los 3 meses de simulación en toda la región.

for execution:
    $ python3 map_mean_precip_wholeperiod.py


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

# rains_acum03_alltimes_wrfout_<dir_of_chunks>.nc
wrfs = [Docpy.WRF(path) for path in glob.glob('./rains*')]
acum = 24                    
lat, lon = wrfs[0].get_latlon()
fail = sys.argv[0][:-3]                                           # name of script
figsdir = os.path.abspath( os.path.join(wrfs[0].ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

### Plot Properties

plt.rcParams.update({'font.size': 14})
box = Docpy.functions.utils['boxHR3']
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='gray',
        facecolor='none',
        )
paleta = plt.get_cmap('gist_ncar')
paleta.set_under('white')
paleta.set_over('magenta')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Calculos
precips = [Docpy.precip.calc_precip_acum(wrf, acum=acum)[:-1] for wrf in wrfs]
# mm/day (tiro el ultimo dia que son 0s)

precips_1mm = [np.where(precip>1, np.ones_like(precip), 0) for precip in precips]
# matriz booleana donde se supera 1mm

frecs_precip1mm = [np.sum(precip_1mm, axis=0)/precip.shape[0] \
        for precip_1mm, precip in zip(precips_1mm, precips)] 
# frecuencia relativa donde se supera la condicion

# Promedio las frecuencias:

frec_precip1mm = np.mean(
        np.stack(
            frecs_precip1mm,
            axis=0,
            ),
        axis=0,
        )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Niveles
levels = np.linspace(0,0.7, 11)

### Plot
fig = plt.figure( figsize=(10,10) )
axes = plt.axes(projection=ccrs.PlateCarree())
mapa = axes.contourf(
        lon, lat, frec_precip1mm,
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

# meto el cbar a mano
fig.subplots_adjust(right=0.8)
cbar_xmin = axes.get_position().xmax+0.05
cbar_ymin = axes.get_position().ymin
cbar_width = 0.025
cbar_hight = axes.get_position().ymax-axes.get_position().ymin
cbar_ax = fig.add_axes([cbar_xmin, cbar_ymin, cbar_width, cbar_hight])
cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
cbar.set_label('Freq')

#plt.show()
sape = '_'.join(
        [
            'freq_precip_above1mmday_OND',
            'acum{:02}'.format(acum),
            wrfs[0].ruta.split('/')[-3],
            wrfs[0].ruta.split('/')[-2],
            ]
        )
fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
plt.close()

