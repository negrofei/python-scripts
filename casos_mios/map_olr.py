"""

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

wrf = Docpy.WRF( sys.argv[1] )
wrf_for_box = Docpy.WRF( sys.argv[2] )
lat, lon = wrf.get_latlon()
fail = sys.argv[0][:-3]                                           # name of script
figsdir = os.path.abspath( os.path.join(wrf.ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

### Plot Properties

plt.rcParams.update({'font.size': 14})
box_plot = [-87.1313, -37.5094, -47.2198, -7.15]
lat_box, lon_box = wrf_for_box.get_latlon()
box_wrfHR = [lon_box.min(), lon_box.max(), lat_box.min(), lat_box.max()]
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='k',
        facecolor='none',
        )
paleta = plt.get_cmap('gist_ncar_r')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Calculos

olr = wrf.get_variables('olr')
label_tiempos,_ = wrf.timedata()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Niveles
levels = np.linspace(100,300, 11)

### Plot
for t in range(olr.shape[0]):
    fig = plt.figure( figsize=(10,10) )
    axes = plt.axes(projection=ccrs.PlateCarree())
    mapa = axes.contourf(
            lon, lat, olr[t,0],
            transform=ccrs.PlateCarree(),
            cmap=paleta,
            extend='both',
            levels=levels,
            )
    axes.set_extent(box_plot)
    axes.add_feature(bordes)
    axes.coastlines('50m', color='k')

    ### Caja
    axes.hlines(y=box_wrfHR[2], xmin=box_wrfHR[0], xmax=box_wrfHR[1], transform=ccrs.PlateCarree())
    axes.hlines(y=box_wrfHR[3], xmin=box_wrfHR[0], xmax=box_wrfHR[1], transform=ccrs.PlateCarree())
    axes.vlines(x=box_wrfHR[0], ymin=box_wrfHR[2], ymax=box_wrfHR[3], transform=ccrs.PlateCarree())
    axes.vlines(x=box_wrfHR[1], ymin=box_wrfHR[2], ymax=box_wrfHR[3], transform=ccrs.PlateCarree())

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

    # Mas configs
    axes.set_title(
            ' '.join(
                [
                    label_tiempos[t],
                    ]
                )
            )

    sape = '_'.join(
            [
                'map_olr',
                wrf.ruta.split('/')[-3],
                wrf.ruta.split('/')[-2],
                wrf.ruta.split('/')[-1],
                wrf.dominio,
                label_tiempos[t],
                ]
            )
    fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
    plt.close()

