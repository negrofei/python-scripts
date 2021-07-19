"""

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
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from tqdm import tqdm
import scipy.ndimage as sci


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Rutas y configs

wrf = Docpy.WRF(sys.argv[1])
wrf_for_box = Docpy.WRF( sys.argv[2] )
lat, lon = wrf.get_latlon()
fail = sys.argv[0][:-3]                                           # name of script
figsdir = os.path.abspath( os.path.join(wrf.ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

### Plot Properties

caso = wrf.ruta.split('/')[-3]
config = wrf.ruta.split('/')[-2]
sim = wrf.ruta.split('/')[-1]

box_plot = [-87.1313, -37.5094, -47.2198, -7.15]
lat_box, lon_box = wrf_for_box.get_latlon()
box_wrfHR = [lon_box.min(), lon_box.max(), lat_box.min(), lat_box.max()]

plt.rcParams.update({'font.size': 14})
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='k',
        facecolor='none',
        )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Calculos

times, delta_t_hs = wrf.timedata() 
every = 3
cada = every//delta_t_hs

# Vientos

plevs = list(wrf.get_variables('lev'))
U = wrf.get_variables('umet')[::cada,plevs.index(850),:,:]
V = wrf.get_variables('vmet')[::cada,plevs.index(850),:,:]

# z850
z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)[::cada,:,:]
smooth = True
if smooth == True:
    a = -8.75
    b = 1.55
    r = (np.diff(lon)[0]+np.diff(lat)[0])/2
    f = lambda x: a*x+b
    z850 = sci.filters.gaussian_filter(z850, f(r))

z850 = np.ma.masked_greater(z850, 2000)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Voy a poner las flechas donde las pone ERA-Interim
intervalo = 2
era = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/ERAI_pl200503_129-131-132-133.nc')
latera, lonera = era.get_latlon()
imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box_plot)
lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])


vlevels = np.linspace(-30,30,9)
zlevels = np.arange(1400,1575,25)
zpaleta = plt.get_cmap('PuBuGn_r')

### Plot
for t in tqdm(range(len(times[::cada])), desc='V + z850'):
    fig = plt.figure(figsize=(10,10))
    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.set_extent(box_plot)
    
    # z850
    mapa = axes.contourf(lon, lat, z850[t],
                        transform=ccrs.PlateCarree(),
                        cmap=zpaleta,
                        extend='both',
                        levels=zlevels)

    contorno = axes.contour(lon, lat, z850[t],
                            transform=ccrs.PlateCarree(),
                            colors='darkgray',
                            levels=zlevels,
                            alpha=0.8,
                            )
    # Vientos
    n = 15
    flechas = axes.quiver(
            lon[lista_lon[::intervalo]], 
            lat[lista_lat[::intervalo]],
            U[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
            V[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
            transform=ccrs.PlateCarree(),
            scale=5, scale_units='inches')
    
    axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')

    axes.add_feature( bordes )
    axes.coastlines('50m', color='k')
    
    # Caja
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
    cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=zlevels)#, format=OOMFormatter(-3))
    cbar.set_label(r'm')
    
    # Mas configs
    axes.set_title(times[::int(every/delta_t_hs)][t], fontsize=14)

    # Savefig
    sape = '_'.join(
            [
                'z850_UV',
                caso,
                config,
                sim,
                wrf.dominio,
                times[::int(every/delta_t_hs)][t],
                ]
            )
    fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
    plt.close()
