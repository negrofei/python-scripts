"""

"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Importo Paquetes

import os
import sys
import pathlib
import numpy as np
import Docpy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import datetime
import matplotlib.dates as mdates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Rutas y configs

ruta = pathlib.Path('.').absolute()
config = 'pblMYNN'
caso = ruta.name
aaaamm = caso[4:8]+caso[9:]

fail = sys.argv[0][:-3]                                           # name of script
figsdir = os.path.abspath( os.path.join(ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

wrf_for_box = Docpy.WRF('/home/martin.feijoo/casos_mios/caso2005-03/pbl_MYNN/4km_B/wrfout_arw_d01_2005-03.nc')

wrf_20km = Docpy.WRF(
        ruta.joinpath(
            '_'.join([
                'wrfout',
                'd01',
                aaaamm,
                config,
                '4A20B',
                'eraigrid.nc'
                ])
            )
        )

wrf_4A20 = Docpy.WRF(
        ruta.joinpath(
            '_'.join([
                'wrfout',
                'd02',
                aaaamm,
                config,
                '4A20B',
                'eraigrid.nc'
                ])
            )
        )

wrf_4B = Docpy.WRF(
        ruta.joinpath(
            '_'.join([
                'wrfout',
                'd01',
                aaaamm,
                config,
                '4B',
                'eraigrid.nc'
                ])
            )
        )

wrf_HR3 = Docpy.WRF(
        ruta.joinpath(
            '_'.join([
                'wrfout',
                'd01',
                aaaamm,
                config,
                'HR3',
                'eraigrid.nc'
                ])
            )
        )
erai = Docpy.ERAI(list(ruta.glob('ERAI*'))[0])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Calculos

# Agarro los tiempos y las lat,lon
times_wrf, delta_t_wrf = wrf_20km.timedata()
every = 6
cada = int(every//delta_t_wrf)
times_era, delta_t_erai = erai.timedata()
times_wrf = times_wrf[::cada]

times = Docpy.functions.common_times([times_era,times_wrf])

idxs_wrfs = times_wrf.index(times[0]), times_wrf.index(times[-1])+1
idxs_era = times_era.index(times[0]), times_era.index(times[-1])+1
lat, lon = erai.get_latlon()

# Agarro el z850 para cada coso

z850_20km = Docpy.campos.calc_z_lvl(wrf_20km, lvl=850)[idxs_wrfs[0]:idxs_wrfs[1],...]
z850_4A20 = Docpy.campos.calc_z_lvl(wrf_4A20, lvl=850)[idxs_wrfs[0]:idxs_wrfs[1],...]
z850_4B = Docpy.campos.calc_z_lvl(wrf_4B, lvl=850)[idxs_wrfs[0]:idxs_wrfs[1],...]
z850_HR3 = Docpy.campos.calc_z_lvl(wrf_HR3, lvl=850)[idxs_wrfs[0]:idxs_wrfs[1],...]
z850_era = Docpy.campos.calc_z_lvl(erai, lvl=850)[idxs_era[0]:idxs_era[1],...]

# Calculo el RSME
# primero obtengo los limites

lat_box, lon_box = wrf_for_box.get_latlon()
imin = int(np.argwhere( np.abs(lon-lon_box[0]) == np.min(np.abs(lon-lon_box[0]))))
imax = int(np.argwhere( np.abs(lon-lon_box[-1]) == np.min(np.abs(lon-lon_box[-1]))))+1
jmin = int(np.argwhere( np.abs(lat-lat_box[0]) == np.min(np.abs(lat-lat_box[0]))))
jmax = int(np.argwhere( np.abs(lat-lat_box[-1]) == np.min(np.abs(lat-lat_box[-1]))))+1

lonbox = lon[imin:imax]
latbox = lat[jmin:jmax]

# corto los campos a la region de HR2

z850_20km = z850_20km[:,jmin:jmax, imin:imax]
z850_4A20 = z850_4A20[:,jmin:jmax, imin:imax]
z850_4B = z850_4B[:,jmin:jmax, imin:imax]
z850_HR3 = z850_HR3[:,jmin:jmax, imin:imax]
z850_era = z850_era[:,jmin:jmax, imin:imax]


# Ahora calculo el rmse para losp untos en el borde y los puntos en el centro

def rmse_borde_centro(field_true, field_pred, ancho_borde):

    # hago el split
    # Para eso me construyo una matriz de 1s en el borde y 0s adentro
    bordes = np.zeros_like(field_true)

    bordes[:,:ancho_borde,:] = 1
    bordes[:,-ancho_borde:,:] = 1
    bordes[:,:,:ancho_borde] = 1
    bordes[:,:,-ancho_borde:] = 1 

    # me armo el inner domain
    inner = -1*bordes+1
    # multiplico las matrices
    bordes_true = np.ma.masked_equal(bordes*field_true, 0)
    inner_true = np.ma.masked_equal(inner*field_true,0)
    
    bordes_pred = np.ma.masked_equal(bordes*field_pred, 0)
    inner_pred = np.ma.masked_equal(inner*field_pred, 0)

    # Calculo el rmse:
    rmse_bordes = [
            np.sum(
                np.sqrt(
                    (bordes_true[t].ravel()-bordes_pred[t].ravel())**2
                    )
                )/(np.sum(bordes[t])) for t in range(len(times))]
    rmse_inner = [
            np.sum(
                np.sqrt(
                    (inner_true[t].ravel()-inner_pred[t].ravel())**2
                    )
                )/(np.sum(inner[t])) for t in range(len(times))]

    return rmse_bordes, rmse_inner

    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot
datetimes = [datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S') for time in times]
plt.rcParams.update({'font.size': 14})
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='k',
        facecolor='none',
        )
box_wrfHR = [lon_box.min(), lon_box.max(), lat_box.min(), lat_box.max()]

for spinup in range(7):
    for ancho in range(1,11):
        print('spinup {0}, ancho {1}'.format(spinup,ancho))

        rmse_borde20, rmse_inner20 = rmse_borde_centro(z850_era,z850_20km,ancho)
        rmse_borde4A20, rmse_inner4A20 = rmse_borde_centro(z850_era,z850_4A20,ancho)
        rmse_borde4B, rmse_inner4B = rmse_borde_centro(z850_era,z850_4B,ancho)
        rmse_bordeHR3, rmse_innerHR3 = rmse_borde_centro(z850_era,z850_HR3,ancho)

        fig, axes = plt.subplots(nrows=1, figsize=(10,10), sharex=True, sharey=True )
        ax2 = axes
        ax2.set_title('Borde de espesor {}'.format(ancho))
        
        ax2.plot(datetimes[spinup:], rmse_inner20[spinup:], label='20kmLR2')
        ax2.plot(datetimes[spinup:], rmse_inner4A20[spinup:], label='4n20')
        ax2.plot(datetimes[spinup:], rmse_inner4B[spinup:], label='4noNest')
        ax2.plot(datetimes[spinup:], rmse_innerHR3[spinup:], label='HR3', color='tab:red')
        
        #ax2.set_title('Resto del dominio')
        for ax in [axes]:
            ax.set_ylabel('RMSE')
            ax.set_yticks(range(5,35,5))
            ax.set_ylim(bottom=2,top=32)
            ax.set_xlim(datetimes[spinup]-datetime.timedelta(hours=3), datetimes[-1]+datetime.timedelta(hours=3))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
            #ax.xaxis.set_minor_locator(mdates.HourLocator(interval=3))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H'))
            ax.set_xticks(datetimes[spinup:])
            ax.grid()

        fig.autofmt_xdate()
        fig.legend(ncol=4, bbox_to_anchor=(0.825,0.96))
        fig.suptitle('ERAI', fontsize=18)


        sape = '_'.join(
                [
                    'rmse',
                    'geopt',
                    'spinup{}'.format(spinup),
                    'espesor{}'.format(ancho),
                    caso,
                    config,
                    'erai',
                    ]
                )
        fig.savefig(os.path.join(figsdir, sape+'.png'), dpi=300, bbxo_inches='tight')
        plt.close()


