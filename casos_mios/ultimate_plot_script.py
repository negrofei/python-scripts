#%% Importo Paquetes

import os
import sys
import numpy as np
import Docpy
import time
##
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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

# props
caso = wrf.ruta.split('/')[4]
config = wrf.ruta.split('/')[5]
sim = wrf.ruta.split('/')[6]
provincias = shpreader.Reader('/home/martin.feijoo/provincias/provincias.shp')
plt.rcParams.update({'font.size': 14})
ruta_figs = os.path.join(wrf.ruta, 'figuras')
try:
    os.mkdir(ruta_figs)
except:
    print('ya esisteeeeee')

while True:
    #opt = input('BOX4 o BOX4B o HR3?? ([1]/2/3) \t')
    opt = '2'
    if (opt.lower() == '1') or (opt == ''):
        box = Docpy.functions.utils['box']
        box_name = 'box4_chica'
        break
    elif opt.lower() == '2':
        box = [-87.1313, -37.5094, -47.2198, -7.15]
        # esta caja es LR2
        box_name = 'box4_grande'
        break
    elif opt.lower() == '3':
        #box = [-76.5693, -47.7057, -36.535, -15.3772]
        box = Docpy.functions.utils['boxHR3']
        # esta caja es HR3
        box_name = 'box_HR3'
        break
    else:
        print('Escribí bien, pelotudo')
#%% Control de figuras
#while True:
#    a = input('Hacemos figuras de precip? ([Y]/n]) \t')
#    if (a.lower()=='y') or (a==''):
#        a='v'
#        break
#    elif a.lower()=='n':
#        break
#    else:
#        print('Escribí bien, pelotudo')
#while True:
#    b = input('Hacemos figuras de Z850? ([Y]/n]) \t')
#    if (b.lower()=='y') or (b==''):
#        b='v'
#        break
#    elif b.lower()=='n':
#        break
#    else:
#        print('Escribí bien, pelotudo')
#
#while True:
#    c = input('Hacemos figuras de Vq + div(Vq)? ([Y]/n]) \t')
#    if (c.lower()=='y') or (c==''):
#        c='v'
#        break
#    elif c.lower()=='n':
#        break
#    else:
#        print('Escribí bien, pelotudo')
#
#while True:
#    d = input('Hacemos figuras de V + |v| + z850)? ([Y]/n]) \t')
#    if (d.lower()=='y') or (d==''):
#        d='v'
#        break
#    elif d.lower()=='n':
#        break
#    else:
#        print('Escribí bien, pelotudo')

#a=b=c=d='v'
a=b=c=d='n'
#e = 'v'
e = 'v'
########################################%% PLOT PRECIP %%############################################
if a=='v':
    paleta = plt.get_cmap('gist_ncar')
    paleta.set_under('white')
    try:
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
                axes.set_extent(box)
                bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                        category='cultural',
                        name='admin_0_boundary_lines_land',
                        scale='50m',
                        edgecolor='k',
                        facecolor='none'
                        )
                axes.add_feature( bordes )
                axes.coastlines('50m')
    #            for rec in provincias.records():
    #                axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
                sape = '_'.join(['acum{:02}'.format(acum), wrf.dominio, caso, config, sim, times[t]])
                fig.savefig(os.path.join('.',sape+'.png'), dpi=300, bbox_inches='tight')
                plt.close()
    except:
        None

########################################%% PLOT Z850 %% ############################################
if b == 'v':
    paleta = plt.get_cmap('cool')#plt.get_cmap('YlGnBu_r')
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
        contorno = axes.contour(lon, lat, z850[t,:,:],
                        transform=ccrs.PlateCarree(),
                        levels=levels)
        axes.set_extent(box)
        bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='50m',
                edgecolor='k',
                facecolor='none'
                )
        axes.add_feature( bordes )
        axes.coastlines('50m')
#        for rec in provincias.records():
#            axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
        sape = '_'.join(['z850', box_name ,wrf.dominio, caso, config, sim, times[t]])
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
    imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box)
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


        axes.set_extent(box)
        bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='50m',
                edgecolor='k',
                facecolor='none'
                )
        axes.add_feature( bordes )
        axes.coastlines('50m')
        #for rec in provincias.records():
        #    axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
        sape = '_'.join(['Vq_cfhi',box_name, wrf.dominio, caso, config, sim, times[::int(every/delta_t_hs)][t]])
        fig.savefig(os.path.join('.',sape+'.png'), dpi=300, bbox_inches='tight')
        plt.close()
        

#############################%% Z850 (cont) + UV (vector) + V (shaded) ##############################
if d == 'v':
    import scipy.ndimage as sci
    times, delta_t_hs = wrf.timedata() 
    every = 3
    cada = every//delta_t_hs
    # Vq
    Vq = Docpy.campos.calc_Vq_integrated(wrf, delta_t=every)
    # vientos 
    plevs = list(wrf.get_variables('lev'))
    #U = wrf.get_variables('umet')[::cada,plevs.index(850),:,:]
    V = wrf.get_variables('vmet')[::cada,plevs.index(850),:,:]
    vpaleta = plt.get_cmap('Spectral_r')
    vlevels = np.linspace(-30,30,9)
    levels = np.arange(1400,1575,25)
    z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)
    # z850
    zlevels = np.arange(1400,1575,25)
    z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)[::cada,:,:]
    smooth = True
    if smooth == True:
        a = -8.75
        b = 1.55
        r = (np.diff(lon)[0]+np.diff(lat)[0])/2
        f = lambda x: a*x+b
        z850 = sci.filters.gaussian_filter(z850, f(r))
    z850 = np.ma.masked_greater(z850, 2000)
    # voy a poner las flechas donde las pone ERA-Interim
    if opt == '2':
        intervalo = 2
    else:
        intervalo = 1
    era = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/UV_ERAI_200503.nc')
    latera, lonera = era.get_latlon()
    imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box)
    lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
    lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])

    # voy a poner algunos contornos rojos para la precipitacion
    # aca tengo 2 opciones: o poner la precipitación acumulada cada 3 horas con un lag
    #                       o poner la precipitacion acumulada cada 6
    # voy a hacer las 2
    levels_pp = [20]
    for opt in ['acum_in_{cada}'.format(cada=cada),'acum_in_running']:
        if opt == 'acum_in_{cada}'.format(cada=cada):
            precip = Docpy.precip.calc_precip_acum(wrf, acum=every)
            lag = 1
        elif opt == 'acum_in_running':
            precip = Docpy.precip.calc_running_acum(wrf, acum=every, run=every*2)
            lag = 0
        for t in tqdm(range(len(times[::cada])-lag-1), desc='V + z850'):
            fig = plt.figure(figsize=(10,10))
            axes = plt.axes(projection=ccrs.PlateCarree())
            axes.set_extent(box)
            mapa = axes.contourf(lon, lat, V[t],
                                transform=ccrs.PlateCarree(),
                                cmap=vpaleta,
                                extend='both',
                                levels=vlevels)
            n = 200
            flechas = axes.quiver(lon[lista_lon[::intervalo]], lat[lista_lat[::intervalo]], 
                    Vq[0][:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                    Vq[1][:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                                  transform=ccrs.PlateCarree(),
                                  scale=5, scale_units='inches')
            axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{kg}{m s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')

            contorno = axes.contour(lon, lat, z850[t],
                                    transform=ccrs.PlateCarree(),
                                    colors='k',
                                    levels=zlevels,
                                    )
            axes.clabel(contorno, inline=True, fmt='%0d')

            contorno2 = axes.contour(lon, lat, precip[t+lag,:,:],
                                    transform=ccrs.PlateCarree(),
                                    colors='purple',
                                    levels=levels_pp,
                                    alpha=0.8,
                                    )
            axes.clabel(contorno2, inline=True, fmt='%0d')
            bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                    category='cultural',
                    name='admin_0_boundary_lines_land',
                    scale='50m',
                    edgecolor='gray',
                    facecolor='none'
                    )
            axes.add_feature( bordes )
            axes.coastlines('50m', color='gray')
            #for rec in provincias.records():
            #    axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
            cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=vlevels)#, format=OOMFormatter(-3))
            cbar.set_label(r'm s$^{-1}$')
            
            #plt.show()
            sape = '_'.join(['Vq_cfhi',opt, box_name, wrf.dominio, caso, config, sim, times[::int(every/delta_t_hs)][t]])
            fig.savefig(os.path.join('.',sape+'.png'), dpi=300, bbox_inches='tight')
            plt.close()

#################### z850 + V + PP20 ######################################################
if e == 'v':
    import scipy.ndimage as sci
    times, delta_t_hs = wrf.timedata() 
    every = 3
    cada = every//delta_t_hs
    # vientos 
    plevs = list(wrf.get_variables('lev'))
    U = wrf.get_variables('umet')[::cada,plevs.index(850),:,:]
    V = wrf.get_variables('vmet')[::cada,plevs.index(850),:,:]
    vpaleta = plt.get_cmap('Spectral_r')
    vlevels = np.linspace(-30,30,9)
    z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)
    # z850
    zlevels = np.arange(1400,1575,25)
    z850 = Docpy.campos.calc_z_lvl(wrf, lvl=850)[::cada,:,:]
    smooth = True
    if smooth == True:
        a = -8.75
        b = 1.55
        r = (np.diff(lon)[0]+np.diff(lat)[0])/2
        f = lambda x: a*x+b
        z850 = sci.filters.gaussian_filter(z850, f(r))
    z850 = np.ma.masked_greater(z850, 2000)
    zpaleta = plt.get_cmap('cool')
    # voy a poner las flechas donde las pone ERA-Interim
    if opt == '2':
        intervalo = 2
    else:
        intervalo = 1
    era = Docpy.ERAI('/home/martin.feijoo/ERA-Interim/casos_mios/ERAI_pl200503_129-131-132-133.nc')
    latera, lonera = era.get_latlon()
    imin, imax, jmin, jmax = Docpy.functions.latlon_idx(latera, lonera, box)
    lista_lat = np.array([np.where(np.abs(lat-j)==np.min(np.abs(lat-j)))[0][0] for j in latera[jmin:jmax+1]])
    lista_lon = np.array([np.where(np.abs(lon-i)==np.min(np.abs(lon-i)))[0][0] for i in lonera[imin:imax+1]])

    # voy a poner algunos contornos rojos para la precipitacion
    # aca tengo 2 opciones: o poner la precipitación acumulada cada 3 horas con un lag
    #                       o poner la precipitacion acumulada cada 6
    # voy a hacer las 2
    #levels_pp = [20]
    #opt = 'acum_in_running'
    #precip = Docpy.precip.calc_running_acum(wrf, acum=every, run=every*2)
    lag = 0
    for t in tqdm(range(len(times[::cada])-lag-1), desc='V + z850'):
        fig = plt.figure(figsize=(10,10))
        axes = plt.axes(projection=ccrs.PlateCarree())
        axes.set_extent(box)
        mapa = axes.contourf(lon, lat, z850[t],
                            transform=ccrs.PlateCarree(),
                            cmap=zpaleta,
                            extend='both',
                            levels=zlevels)
        n = 15
        flechas = axes.quiver(lon[lista_lon[::intervalo]], lat[lista_lat[::intervalo]], 
                U[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                V[:,lista_lat[::intervalo],:][:,:,lista_lon[::intervalo]][t]/n,
                              transform=ccrs.PlateCarree(),
                              scale=5, scale_units='inches')
        axes.quiverkey(flechas,.875,1.025, 1, '{}'.format(n)+r'$\frac{m}{s}$', labelpos='E', coordinates='axes', fontproperties={'size':15}, color='k', labelcolor='k')

        contorno = axes.contour(lon, lat, z850[t],
                                transform=ccrs.PlateCarree(),
                                colors='darkgray',
                                levels=zlevels,
                                alpha=0.8,
                                )
        #axes.clabel(contorno, inline=True, fmt='%0d')

        #contorno2 = axes.contour(lon, lat, precip[t+lag,:,:],
        #                        transform=ccrs.PlateCarree(),
        #                        colors='red',
        #                        levels=levels_pp,
        #                        alpha=0.8,
        #                        )
        #axes.clabel(contorno2, inline=True, fmt='%0d')
        bordes = cfeature.NaturalEarthFeature(          # limites de los paises
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='50m',
                edgecolor='k',
                facecolor='none'
                )
        axes.add_feature( bordes )
        axes.coastlines('50m', color='k')
        #for rec in provincias.records():
        #    axes.add_geometries( [rec.geometry], ccrs.PlateCarree(), edgecolor="k", facecolor='none', linewidth=0.5)
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
        cbar = fig.colorbar(mapa, cax=cbar_ax, ticks=zlevels)#, format=OOMFormatter(-3))
        cbar.set_label(r'm')
        #if t == 0 :
            #plt.show()
        sape = '_'.join(['z850_UV',opt, box_name, wrf.dominio, caso, config, sim, times[::int(every/delta_t_hs)][t]])
        fig.savefig(os.path.join(ruta_figs,sape+'.png'), dpi=300, bbox_inches='tight')
        plt.close()
