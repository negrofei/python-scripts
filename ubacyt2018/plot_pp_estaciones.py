import Docpy
import sys
import os
import numpy as np
import pandas as pd
#
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ruta = sys.argv[1]
arch = sys.argv[0]
data = Docpy.Dato(ruta)
acum = 24
lat, lon = data.get_latlon()
precip = Docpy.precip.calc_precip_acum(data, acum=acum)

### Agarro las estaciones que se usaron para generar el interpolado con lo que me paso mati de ubacyt

data_stations = pd.read_csv('/home/martin.feijoo/ubacyt2018/dominios_wrf_ESD/StationData.txt',sep=' ').to_numpy()[:,1:]   # agarre las latlon de las estaciones
box = Docpy.functions.utils['box']
#grafico
times = Docpy.functions.acum_time(data, acum=acum)
#provincias = shpreader.Reader('/home/martin.feijoo/provincias/provincias.shp')
plt.rcParams.update({'font.size': 14})
paleta = plt.get_cmap('gist_ncar')
paleta.set_under('white')
#levels = np.arange(35,385, 35)
levels = Docpy.precip.d_format['levels'][acum]
Docpy.functions.printer('Grafico precipitation...')
for t in range(len(times)):
    fig = plt.figure(figsize=(10,10))
    axes = plt.axes(projection=ccrs.PlateCarree())
    mapa = axes.contourf(lon, lat, precip[t,:,:],
                         transform=ccrs.PlateCarree(),
                         cmap=paleta,
                         extend='both',
                         levels=levels)
    axes.set_extent(box)
#    for rec in provincias.records():
#        axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
    
    # Agrego las estaciones
    for station in data_stations:
        if (box[0]<station[1]<box[1]) and (box[2]<station[0]<box[3]):
            axes.scatter(station[1],station[0], color='k', transform=ccrs.PlateCarree(), s=2)

    gl = axes.gridlines(draw_labels=True, color='black', alpha=0.2, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':14, 'rotation':25}
    gl.ylabel_style = {'size':14}
    axes.set_title('estaciones', fontsize=14)

    # meto un cbar a mano
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
    cbar.set_label('mm')
    if t == 0:
        plt.show()
    ##
    Docpy.functions.printer('SAVEFIG---')
    sape = '_'.join(['precip',data.filename.split('.')[0],'acum{:02d}'.format(acum), times[t][:-9]])
    fig.savefig(os.path.join('.', sape+'.jpg'), dpi=150, bbox_inches='tight')
    plt.close()
Docpy.functions.printer('LISTORTI')

