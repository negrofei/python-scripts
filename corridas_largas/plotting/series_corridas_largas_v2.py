"""
Grafico series temporales de la precipitación de una determinada corrida larga junto con las 3 bases
de datos que tenemos siempre: cmorph, imerg y mswep

for execution:
    $ python3 series_corridas_largas_v1.py 

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
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from tqdm import tqdm
import re

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
### Rutas y configs

wrfs = [Docpy.WRF(path) for path in glob.glob('./rains*')]
# rains_acum03_alltimes_wrfout_<dir_of_chunks>.nc

observaciones = [               # agarro los datos de 2015
        Docpy.Dato(glob.glob('/home/martin.feijoo/CMORPH/casos_mios/200501-201612_8km/*2015*')[0]),
        Docpy.Dato(glob.glob('/home/martin.feijoo/GPM/*2015*')[0]),
        Docpy.Dato(glob.glob('/home/martin.feijoo/MSWEP/0.1deg_3hly_data/*2015*')[0]),
        ]
print('Se graficaran las series de lo siguiente',         # informo que archivos agarro
        os.path.join(wrfs[0].ruta,wrfs[0].filename)+'\n', sep='\n')
[print(os.path.join(obs.ruta,obs.filename)+'\n') for obs in observaciones]

acum_list = [24, 6, 3]          # grafico para acumular cada 3 6 y 24
lat, lon = wrfs[0].get_latlon()
fail = sys.argv[0][:-3]         # name of scripts 
figsdir = os.path.abspath( os.path.join(wrfs[0].ruta,'figuras',fail) )
try:
    os.makedirs( figsdir )
except:
    None

### Plot Properties

plt.rcParams.update({'font.size': 14})
paleta = plt.get_cmap('gist_ncar')
paleta.set_under('white')
linestyles_obs = ['-','--','-.']
acum_params = {
        24:{'yticks':np.arange(0,45,5),
            'x_every':3},
        6:{'yticks':np.arange(0,20,2.5),
            'x_every':12},
        3:{'yticks':np.arange(0,14,2),
            'x_every':24},
        }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Mean precipitation calculation ###
for acum in acum_list:                  # hago para todas las acumulaciones
    Docpy.functions.printer('acum{:02}'.format(acum))
    box = Docpy.functions.utils['box']  # DATO DE VITAL IMPORTANCIA: agarro esta region para la caja

    # OBS
    ppmean_obs = []
    time_obs = []
    for obs in observaciones:
        Docpy.functions.printer('Mean PP',obs.filename)
        ppmean_obs.append( np.nan_to_num( Docpy.precip.ppmean_box(obs, acum, box, offset=0) ) )
        time_obs.append( Docpy.functions.acum_time(obs, acum) )

    # WRF
    ppmean_wrf = []
    time_wrf = []
    for worf in wrfs:
        Docpy.functions.printer('Mean PP', worf.filename)
        ppmean_wrf.append( np.nan_to_num( Docpy.precip.ppmean_box(worf, acum, box, offset=0) ) )
        time_wrf.append( worf.acum_time(acum) )

    # Common times between OBS and WRF
    time_common = Docpy.functions.common_times(time_obs + time_wrf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Plot

    fig, ax = plt.subplots(figsize=(15,7))
    #WRF
    # Ensemble
    sliwrf_ens = slice(time_wrf[0].index(time_common[0]), time_wrf[0].index(time_common[-1]))
    mean = np.mean(
            np.stack(
                ppmean_wrf,
                axis=0,
                ),
            axis=0,
            )
    ax.plot(time_wrf[0][sliwrf_ens], mean[sliwrf_ens],
            color = 'r',
            linestyle='solid',
            linewidth=1,
            label=wrfs[0].name)
    # Members
    for w, worf in enumerate(wrfs):
        sliwrf = slice(time_wrf[w].index(time_common[0]), time_wrf[w].index(time_common[-1]))
        ax.plot(time_wrf[w][sliwrf], ppmean_wrf[w][sliwrf],
                color='r',
                linestyle='dashed',
                linewidth=1,
                )


    #OBS
    # Correlation coeficient
    rs = []
    for o, obs in enumerate(observaciones):
        sliobs = slice(time_obs[o].index(time_common[0]), time_obs[o].index(time_common[-1]))
        
        rs.append(np.corrcoef(x=ppmean_wrf[0][sliwrf_ens], y=ppmean_obs[o][sliobs])[0,1])
    
        ax.plot(time_obs[o][sliobs], ppmean_obs[o][sliobs],
                color='k',
                linestyle=linestyles_obs[o],
                linewidth=1,
                label=obs.name+' ('+str(round(rs[o],3))+')')
        # Plot single events
        if acum == 24:
            idx = time_obs[o].index('2015-11-10 00:00:00')
            ax.scatter(time_obs[o][idx], ppmean_obs[o][idx], color='blue')

    # Config y
    yticks = acum_params[acum]['yticks']          # acomodado a mano
    ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
    ax.set_yticks(list(yticks))
    ax.set_yticklabels(
            ['{}'.format(tick) for tick in ax.get_yticks()]
            )
    ax.set_ylabel('mm')
    # Config x
    x_every = acum_params[acum]['x_every']
    ax.set_xticks(time_common[::x_every])
    ax.set_xticklabels([time[:-6] for time in time_common][::x_every], fontsize=14)
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    
    # Mas configs
    ax.grid(alpha=0.7)
    ax.legend()
    ax.set_title(' '.join(['Mean precipitation',
        os.path.split(wrfs[0].ruta)[-1],
        'box4',
        'acum{:02d}'.format(acum),
        ]
        )
        )

    # Savefig
    sape = '_'.join(
            [
                'meanpp',
                'box4',
                'acum{:02}'.format(acum),
                os.path.split(wrfs[0].ruta)[-1],
                ]
            )
    fig.savefig(os.path.join(figsdir,sape+'.png'), dpi=300, bbox_inches='tight')
    plt.close()

