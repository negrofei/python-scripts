#!/home/nadia/anaconda3/bin/python3

"""
Script para calcular la precipitación media en una caja comun entre simulaciones y observaciones
"""

import os
import sys

import Docpy
from Docpy import Dato
from Docpy import WRF
import matplotlib.pyplot as plt

import numpy as np
###### cosas a mano

fail = sys.argv[0].split('.')[0]

#### rutas

ruta = os.path.abspath('./')
ruta_figs = os.path.join(ruta, 'figuras')
try: 
    os.mkdir(ruta_figs)
    print(ruta_gis, 'creado con exito')
except:
    print(ruta_figs, 'ya existe')

# definiciones
res_linestyle = {'2.4km':'dashed',
                 '4km':'solid',
                 '12km':'dashdot',
                 '20km':'dotted'}

simus = ['2.4km_anidado12',
         '4km_anidado12',
         '4km_anidado12_B',
         '4km_anidado20_B',
         '4km_B',
         '4km_B_15S47W']

configs = [
        'config_inicial',
        'config_inicialYSU',
        'igual1_con-shaYSU',
        'shallow-on_MYNN_grell_MP',
        ]

l = len(configs)

domains = dict(zip(simus,['1A12','1A12','2B12','2B20','2B','3']))

size_colors = {'1A12':'tab:blue',
               '2B12':'tab:orange',
               '2B20':'chocolate',
               '2B'  :'goldenrod',
               '3':'tab:red'}
size_colors_names = dict(zip([*size_colors], ['1-n12','2-n12','2-n20','2-noNest','3-noNest']))

caso = ruta.split('/')[['caso2' in x for x in ruta.split('/')].index(True)]
datos = ['cmorph','mswep','imerg']
casito = caso[4:]


####
wrfs = {}
for sim in simus:
    if sim==[*simus][-1]:
        dom = 'd01'
        for config in configs:
            ruta_wrf = os.path.abspath(
                    os.path.join(
                        ruta,
                        '..',
                        '..',
                        config,
                        sim,
                        '_'.join(['wrfout_arw_d01',casito+'.nc'])
                        )
                    )
            try:
                wrf = WRF(ruta_wrf)
            except:
                print(ruta_wrf,'not found')
            else:
                res = Docpy.functions.get_namelist_info(os.path.join(wrf.ruta,'namelist.input'))[0]
                wrf.color = size_colors[domains[sim]]
                r = res[int(dom[-1])-1]
                if r == 2400:
                    wrf.linestyle = res_linestyle[str(r/1000)+'km']
                else:
                    wrf.linestyle = res_linestyle[str(int(r/1000))+'km']

                if r  > 4000:
                    wrf.size = 'LR'
                else:
                    wrf.size = 'HR'
                wrf.sim = sim
                if r == 2400:
                    wrf.label = str(r/1000)+'km'+'_'+wrf.size+size_colors_names[domains[sim]]
                else:
                    wrf.label = str(int(r/1000))+'km'+'_'+wrf.size+size_colors_names[domains[sim]]
                wrfs[config,sim,dom] = wrf
    config = ruta.split('/')[5]
    for dom in ['d01','d02']:
        ruta_wrf = os.path.abspath(os.path.join(ruta,'..',sim,'_'.join(['wrfout',
                                                                        'arw',
                                                                        dom,
                                                                        casito+'.nc'])))
        try: 
            wrf = WRF(ruta_wrf)
        except:
            print(ruta_wrf,'not found')
        else:
            res = Docpy.functions.get_namelist_info(os.path.join(wrf.ruta,'namelist.input'))[0]
            wrf.color = size_colors[domains[sim]]                   
            r = res[int(dom[-1])-1]
            if r == 2400:
                wrf.linestyle = res_linestyle[str(r/1000)+'km']
            else:
                wrf.linestyle = res_linestyle[str(int(r/1000))+'km']

            if r  > 4000:
                wrf.size = 'LR'
            else:
                wrf.size = 'HR'
            wrf.sim = sim
            if r == 2400:
                wrf.label = str(r/1000)+'km'+'_'+wrf.size+size_colors_names[domains[sim]]
            else:
                wrf.label = str(int(r/1000))+'km'+'_'+wrf.size+size_colors_names[domains[sim]]
            wrfs[config, sim, dom] = wrf
#sys.exit()
#### chequeo que esto haya andado
for wrf in [*wrfs]:
    print(wrfs[wrf].sim,wrfs[wrf].label)

observaciones = [
        Dato('/home/martin.feijoo/CMORPH/casos_mios/'+caso[4:8]+caso[9:11]+'/CMORPH_V1.0_ADJ_8km-30min_'+caso[4:8]+caso[9:11]+'_CSAM1.nc'),
        Dato('/home/martin.feijoo/MSWEP/0.1deg_3hly_data/'+caso[4:8]+caso[9:11]+'_CSAM.nc'),
        Dato('/home/martin.feijoo/GPM/'+caso[4:8]+caso[9:11]+'/3B-HHR.MS.MRG.3IMERG.'+caso[4:8]+caso[9:11]+'.nc'),
        ]
##### 

acum = int(input('Acumulacion de precipitacion? '))
box = Docpy.functions.common_doms([os.path.join(wrfs[wrf].ruta,wrfs[wrf].filename) for wrf in [*wrfs]])

# tiro los archivos que no cumplen el DT de acumulacion
for obs in observaciones:
    if acum < obs.timedata()[1]:
        observaciones.pop(observaciones.index(obs))

import sys
sys.exit()
# pp mean
ppmean_obs = []
time_obs = []
for obs in observaciones:
    ppmean_obs.append( Docpy.precip.ppmean_box(obs, acum, box, offset=.0) )  # UPDATE 13/10/20: promedio los puntos que superan el 1mm / acum hs 
    time_obs.append( Docpy.functions.acum_time(obs, acum) )
ppmean_wrf = {}
time_wrf = []
for wrf in [*wrfs]:
    ppmean_wrf[wrf] = Docpy.precip.ppmean_box(wrfs[wrf], acum, box, offset=.0)   # UPDATE 13/10/20: promedio los puntos que superan el 1mm / acum hs
    time_wrf.append( wrfs[wrf].acum_time(acum) )

# me agarro los tiempos comunes
time_common = Docpy.functions.common_times(time_obs+time_wrf)

####### GRAFICO
linestyles = ['-', '--', '-.', ':']
plt.rcParams.update({'font.size':14})
fig, ax = plt.subplots(figsize=(10,7))
for o, obs in enumerate(observaciones):
    sli = slice(time_obs[o].index(time_common[0]), time_obs[o].index(time_common[-1]))
    ax.plot(time_obs[o][sli], ppmean_obs[o][sli],
            color='black',
            linestyle=linestyles[o],
            linewidth=2.,
            label=obs.name)
for s, wrf in enumerate([*wrfs][:-l]):
    sli = slice(time_wrf[s].index(time_common[0]), time_wrf[s].index(time_common[-1]))
    ax.plot(time_wrf[s][sli], ppmean_wrf[wrf][sli],
            color=wrfs[wrf].color,
            linestyle=wrfs[wrf].linestyle,
            linewidth=2.,
            label=wrfs[wrf].label)

#for s, wrf in enumerate([*wrfs][-3:]):
#    sli = slice(time_wrf[s].index(time_common[0]), time_wrf[s].index(time_common[-1]))
#    ax.plot(time_wrf[s][sli], ppmean_wrf[wrf][sli],
#            color=k,
#            linestyle=wrfs[wrf].linestyle,
#            linewidth=2.,
#            label=wrfs[wrf].label)   
    
# grafico el fill between esquemas
#sli = slice(time_wrf[s+1].index(time_common[0]), time_wrf[s+1].index(time_common[-1]))

def minimo_de_varios(lista):
    assert len(lista)>=1
    if len(lista)<2:
        return np.minimum(lista[0],lista[0])
    else:
        x = np.minimum(lista[0], lista[1])
        for l in range(2, len(lista)):
            y = np.minimum(x, lista[l])
            x = y
        return x

def maximo_de_varios(lista):
    assert len(lista)>=1
    if len(lista)<2:
        return np.maximum(lista[0],lista[0])
    else:
        x = np.maximum(lista[0], lista[1])
        for l in range(2, len(lista)):
            y = np.maximum(x, lista[l])
            x = y
        return x
lower = minimo_de_varios( [ppmean_wrf[ [*wrfs][-j] ] for j in reversed(range(l))] )
upper = maximo_de_varios( [ppmean_wrf[ [*wrfs][-j] ] for j in reversed(range(l))] )

ax.fill_between(time_wrf[s+1][sli], lower[sli] , upper[sli],
        color=wrfs[[*wrfs][-l]].color,
        linestyle=wrfs[[*wrfs][-l]].linestyle,
        linewidth=2.,
        label=wrfs[[*wrfs][-l]].label,
        )

yticks = np.arange(0,14,2)
ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
ax.set_yticks(list(yticks))
ax.set_yticklabels(
        ['{}'.format(tick) for tick in ax.get_yticks()])
ax.set_ylabel('mm')
ax.set_xticklabels([time_common[i][:-6] for i in range(len(time_common))], fontsize=14)
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
if casito == '2005-03':
    ax.legend(fontsize=14)
ax.grid()
plt.title(casito, fontsize=14)
plt.show()
#%% SAVE FIG
for fmt in ['.jpg']:
    fig.savefig(os.path.join(ruta_figs,
                             '_'.join([fail,caso,'acum{:02d}'.format(acum)]))+fmt, dpi=150, bbox_inches='tight')



