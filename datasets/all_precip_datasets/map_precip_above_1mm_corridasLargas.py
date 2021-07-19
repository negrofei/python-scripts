"""
Mapa de la frecuencia a la cual la precipitación excede el 1mm/day para las bases de datos

for execution:
    $ python3 map_precip_above_1mm_corridasLargas.py


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

data = Docpy.Dato( sys.argv[1] )
acum = 24                    
lat, lon = data.get_latlon()
fail = sys.argv[0][:-3]                       # name of script
figsdir = os.path.abspath( os.path.join(data.ruta,'figuras',fail) )
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
precip = Docpy.precip.calc_precip_acum(data, acum=acum)[:90*int(24/acum)]
# mm/day (tiro el ultimo dia que son 0s)

precip_1mm = np.ma.where(precip>1, np.ones_like(precip), 0) 
# matriz booleana donde se supera 1mm

frec_precip1mm = np.ma.sum(precip_1mm, axis=0)/precip.shape[0] 
# frecuencia relativa donde se supera la condicion


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

# Seteo title
axes.set_title(
        ' '.join(
            [
                'frec precip>1mm/day',
                data.name,
                ]
            )
        )

#plt.show()
sape = '_'.join(
        [
            'freq_precip_above1mmday_OND',
            'acum{:02}'.format(acum),
            data.name,
            ]
        )
fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
plt.close()

