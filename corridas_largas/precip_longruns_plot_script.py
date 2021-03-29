#%% Importo Paquetes

import os
import sys
import numpy as np
import Docpy
import time
##
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
##
from tqdm import tqdm
#%% Rutas y configs
wrf = Docpy.WRF(sys.argv[1])
acum_list = [3, 6, 24] # grafico para acumular cada 3 6 y 24
lat, lon = wrf.get_latlon()
fail = sys.argv[0][:-3]
figsdir = os.path.abspath( os.path.join('.','figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None
# props
provincias = shpreader.Reader('/home/martin.feijoo/provincias/provincias.shp')
plt.rcParams.update({'font.size': 14})


#%% Control de figuras
while True:
    a = input('Hacemos figuras de precip? ([Y]/n]) \t')
    if (a.lower()=='y') or (a==''):
        a='v'
        break
    elif a.lower()=='n':
        break
    else:
        print('Escribí bien, pelotudo')
while True:
    b = input('Hacemos figuras de Z850? ([Y]/n]) \t')
    if (b.lower()=='y') or (b==''):
        b='v'
        break
    elif b.lower()=='n':
        break
    else:
        print('Escribí bien, pelotudo')

while True:
    c = input('Hacemos figuras de Vq + div(Vq)? ([Y]/n]) \t')
    if (c.lower()=='y') or (c==''):
        c='v'
        break
    elif c.lower()=='n':
        break
    else:
        print('Escribí bien, pelotudo')

########################################%% PLOT PRECIP %%############################################
if a=='v':
    paleta = plt.get_cmap('gist_ncar')
    paleta.set_under('white')
    for acum in acum_list:
        levels = Docpy.precip.d_format['levels'][acum]
        precip = Docpy.precip.calc_precip_acum(wrf, acum=acum)
        times = Docpy.functions.acum_time(wrf, acum=acum)
        for t in tqdm(range(len(times)), desc='{:02}'.format(acum)):
            fig = plt.figure(figsize=(10,10))
            axes = plt.axes(projection=ccrs.PlateCarree())
            mapa = axes.contourf(lon, lat, precip[t,:,:],
                            transform=ccrs.PlateCarree(),
                            cmap=paleta,
                            extend='both',
                            levels=levels)
            axes.set_extent(Docpy.functions.utils['box'])
            for rec in provincias.records():
                axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
            gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 14, 'rotation': 25}
            gl.ylabel_style = {'size': 14, 'rotation': 25}

            axes.set_title(times[t], fontsize=14)

            # meto el cbar a mano
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
            cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
            cbar.set_label('mm')

            #plt.show()
            sape = '_'.join(['acum{:02}'.format(acum), wrf.dominio, wrf.ruta.split('/')[2], wrf.ruta.split('/')[3], times[t]])
            fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
            plt.close()

########################################%% PLOT Z850 %% ############################################
if b == 'v':
    paleta = plt.get_cmap('YlGnBu_r')
    levels = np.arange(1400,1575,25)
    z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)
    times, delta_t_hs = wrf.timedata()
    every = 3
    for t in tqdm(range(0,len(times),int(every/delta_t_hs)), desc='Z850'):
        fig = plt.figure(figsize=(10,10))
        axes = plt.axes(projection=ccrs.PlateCarree())
        mapa = axes.contourf(lon, lat, z850[t,:,:],
                        transform=ccrs.PlateCarree(),
                        cmap=paleta,
                        extend='both',
                        levels=levels)
        axes.set_extent(Docpy.functions.utils['box'])
        for rec in provincias.records():
            axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
        gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14, 'rotation': 25}
        gl.ylabel_style = {'size': 14, 'rotation': 25}
        
        axes.set_title(times[t], fontsize=14)
        
        # meto el cbar a mano
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels)
        cbar.set_label('m')
        
        #plt.show()
        sape = '_'.join(['z850', wrf.dominio, caso, config, sim, times[t]])
        fig.savefig(os.path.join('.',sape+'.png'), dpi=300, bbox_inches='tight')
        plt.close()
        

########################################%% Z850 (cont) + Vq (vect) + Conv(Vq) (shaded) #############
if c == 'v':
    import matplotlib.ticker
    class OOMFormatter(matplotlib.ticker.ScalarFormatter):
        def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
            self.oom = order
            self.fformat = fformat
            matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
        def _set_orderOfMagnitude(self, nothing):
            self.orderOfMagnitude = self.oom
            
        def _set_format(self, vmin, vmax):
            self.format = self.fformat
            if self._useMathText:
                self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)
    times, delta_t_hs = wrf.timedata() 
    every = 3
    Vq = Docpy.campos.calc_Vq_integrated(wrf, delta_t=every)
    chif = Docpy.campos.calc_cfhi(wrf, delta_t=every)
    paleta = plt.get_cmap('Spectral_r')
    levels = np.linspace(-1e-2,1e-2,9)
    # voy a poner las flechas donde las pone ERA-Interim
    era = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/UV_ERAI_200503.nc')
    latera, lonera = era.get_latlon()
    imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, Docpy.functions.utils['box'])
    lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
    lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])

    for t in tqdm(range(len(times[::int(every/delta_t_hs)])), desc='Vq+cfhi'):
        fig = plt.figure(figsize=(10,10))
        axes = plt.axes(projection=ccrs.PlateCarree())
        mapa = axes.contourf(lon, lat, chif[t],
                            transform=ccrs.PlateCarree(),
                            cmap=paleta,
                            extend='both',
                            levels=levels)
        n = 200
        flechas = axes.quiver(lon[lista_lon], lat[lista_lat], 
                              Vq[0][:,lista_lat,:][:,:,lista_lon][t]/n,
                              Vq[1][:,lista_lat,:][:,:,lista_lon][t]/n,
                              transform=ccrs.PlateCarree(),
                              scale=5, scale_units='inches')
        axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{kg}{m s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')


        axes.set_extent(Docpy.functions.utils['box'])
        for rec in provincias.records():
            axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
        gl = axes.gridlines(draw_labels=True, color='k', alpha=0.2, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14, 'rotation': 25}
        gl.ylabel_style = {'size': 14, 'rotation': 25}
        
        axes.set_title(times[::int(every/delta_t_hs)][t], fontsize=14)
        
        # meto el cbar a mano
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=levels, format=OOMFormatter(-3))
        cbar.set_label(r'kg s$^{-1}$')
        
        #plt.show()
        sape = '_'.join(['Vq_cfhi', wrf.dominio, caso, config, sim, times[::int(every/delta_t_hs)][t]])
        fig.savefig(os.path.join('.',sape+'.png'), dpi=300, bbox_inches='tight')
        plt.close()
        



