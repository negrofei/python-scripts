#!/home/nadia/anaconda3/bin/python3

"""
Calculo la serie de la precipitacion en la caja comun en los dominios para el ensemble de fisicas.
Es decir, para cada configuracion de dominio ('sim') calculo la envolvente de las físicas (configs)
"""

import os
import sys
import datetime
import Docpy
from Docpy import Dato
from Docpy import WRF
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
###### cosas a mano

fail = sys.argv[0].split('.')[0]

#### rutas

ruta = os.path.abspath('../')
ruta_figs = os.path.join(ruta, 'figuras')
try: 
    os.mkdir(ruta_figs)
    print(ruta_figs, 'creado con exito')
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
        'pbl_MYNN',
        'shallow-on_MYNN_grell',
        'shallow-on_MYNN_grell_MP',
        ]

simus_d0s = dict(zip(simus, [['d02'],['d01','d02'],['d01','d02'],['d01','d02'],['d01'],['d01']]))

l = len(configs)

domains = dict(zip(simus,['1A12','1A12','2B12','2B20','2B','3']))

size_colors = {'1A12':'tab:blue',
               '2B12':'tab:orange',
               '2B20':'chocolate',
               '2B'  :'goldenrod',
               '3':'tab:red'}
size_colors_names = dict(zip([*size_colors], ['1(12)','2(12)','2(20)','2','3']))

caso = ruta.split('/')[['caso2' in x for x in ruta.split('/')].index(True)]
datos = ['cmorph','mswep','imerg']
casito = caso[4:]
name_casito = {
        '2005-03':'Case1 2005-03 ',
        '2015-11':'Case2 2015-11',
        '2016-10':'Case3 2016-10',
        }


####
wrfs = {}
for sim in simus:
    for dom in simus_d0s[sim]:
        for config in configs:
            ruta_wrf = os.path.abspath(
                    os.path.join(
                        ruta,
                        '..',
                        config,
                        sim,
                        '_'.join(['wrfout','arw',dom,casito+'_boxed.nc'])
                        )
                    )
            try:
                wrf = WRF(ruta_wrf)
            except FileNotFoundError as err:
                print(err.strerror,ruta_wrf)
            else:
                wrf.config = config
                # saco el namelist de la carpeta oroginal.
                orig_path = os.path.abspath(wrf.ruta.split('interpolated_files')[0]+wrf.ruta.split('interpolated_files')[1])
                res = Docpy.functions.get_namelist_info(os.path.join(orig_path,'namelist.input'))[0]
                wrf.color = size_colors[domains[sim]] 
                r = res[int(dom[-1])-1]
                if r == 2400:
                    wrf.linestyle = res_linestyle[str(r/1000)+'km']
                else:
                    wrf.linestyle = res_linestyle[str(int(r/1000))+'km']

                if r  > 4000:
                    wrf.size = 'DR'
                else:
                    wrf.size = 'CP'
                wrf.sim = sim
                if r == 2400:
                    wrf.label = wrf.size+size_colors_names[domains[sim]][0]+'-'+str(r/1000)+size_colors_names[domains[sim]][1:]
                else:
                    wrf.label = wrf.size+size_colors_names[domains[sim]][0]+'-'+str(int(r/1000))+size_colors_names[domains[sim]][1:]
                if wrf.size == 'DR':
                    wrf.label = wrf.label[:-4]
                wrfs[config, sim, dom] = wrf
#sys.exit()
#### chequeo que esto haya andado
#for wrf in [*wrfs]:
#    print(wrfs[wrf].config,wrfs[wrf].sim,wrfs[wrf].label)

observaciones = [
        Dato('/home/martin.feijoo/casos_mios/interpolated_files/CMORPH_V1.0_ADJ_8km-30min_'+caso[4:8]+caso[9:11]+'_CSAM1.nc'),
        Dato('/home/martin.feijoo/casos_mios/interpolated_files/mswep_'+caso[4:8]+caso[9:11]+'_boxed.nc'),
        Dato('/home/martin.feijoo/casos_mios/interpolated_files/3B-HHR.MS.MRG.3IMERG.'+caso[4:8]+caso[9:11]+'_boxed.nc'),
        ]
##### 

acum = int(input('Acumulacion de precipitacion? '))
box = Docpy.functions.common_doms([os.path.join(wrfs[wrf].ruta,wrfs[wrf].filename) for wrf in [*wrfs]])

# tiro los archivos que no cumplen el DT de acumulacion
for obs in observaciones:
    if acum < obs.timedata()[1]:
        observaciones.pop(observaciones.index(obs))

# pp mean
ppmean_obs = []
time_obs = []
for obs in observaciones:
    ppmean_obs.append( Docpy.precip.ppmean_box(obs, acum, box, offset=.0) )  # UPDATE 13/10/20: promedio los puntos que superan el offset / acum hs 
    time_obs.append( Docpy.functions.acum_time(obs, acum) )
ppmean_wrf = {}
time_wrf = []
for wrf in [*wrfs]:
    ppmean_wrf[wrf] = Docpy.precip.ppmean_box(wrfs[wrf], acum, box, offset=.0)   # UPDATE 13/10/20: promedio los puntos que superan el offset / acum hs
    # como hay algunos valores altisimos falopa en los últimos tiempos los hago nan
    ppmean_wrf[wrf][-1] = np.nan
    time_wrf.append( wrfs[wrf].acum_time(acum) )


# me agarro los tiempos comunes
time_common = Docpy.functions.common_times(time_obs+time_wrf)

# armo datetimes

time_obs_dt = [[datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S') for time in time_obss] for time_obss in time_obs]

time_wrf_dt = [[datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S') for time in time_wrfss] for time_wrfss in time_wrf]

time_common_dt = [datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S') for time in time_common]


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

####### GRAFICO
linestyles = ['-', '--', '-.', ':']
plt.rcParams.update({'font.size':14})
#for fig_type in ['todasJuntas','split','todasSeparadas','splitNestCP','3WaySplit']:
for fig_type in ['3WaySplit']:
    if fig_type == 'todasJuntas':
        fig, ax = plt.subplots(figsize=(10,7))
        for o, obs in enumerate(observaciones):
            sli = slice(time_obs[o].index(time_common[2]), time_obs[o].index(time_common[-1]))
            ax.plot(time_obs[o][sli], ppmean_obs[o][sli],
                    color='black',
                    linestyle=linestyles[o],
                    linewidth=2.,
                    label=obs.name)

        # grafico el fill between esquemas
        for c in range(round(len([*wrfs])/len(configs))):
            estas_sims = [*wrfs][len(configs)*c:len(configs)*(c+1)]
            sim_name = estas_sims[0][1]
            lower = minimo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
            upper = maximo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
            all_array = np.zeros( (len(estas_sims), len(time_common)) )
            for s, sim in enumerate(estas_sims):
                all_array[s,:] = ppmean_wrf[sim]
            media = np.mean(all_array, axis=0)
            ax.fill_between(time_wrf[0][2:], lower [2:], upper[2:],
                    color=wrfs[estas_sims[0]].color,
                    linestyle=wrfs[estas_sims[0]].linestyle,
                    linewidth=2,
                    #label=wrfs[estas_sims[0]].label,
                    alpha=0.4,
                    )
            ax.plot(time_wrf[0][2:], media[2:],
                    color=wrfs[estas_sims[0]].color,
                    linestyle=wrfs[estas_sims[0]].linestyle,
                    linewidth=2,
                    label=wrfs[estas_sims[0]].label,
                    )

        yticks = np.arange(0,14,2)
        ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
        ax.set_yticks(list(yticks))
        ax.set_yticklabels(
                ['{}'.format(tick) for tick in ax.get_yticks()])
        ax.set_ylabel('mm')
        ax.set_xticklabels([time[:-6] for time in time_common[2:]], fontsize=14)
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        #if casito == '2005-03':
        #    ax.legend(fontsize=14)
        ax.grid()
        plt.title(casito, fontsize=14)
        #plt.show()
        #%% SAVE FIG
        for fmt in ['.jpg']:
            fig.savefig(os.path.join(ruta_figs,
                                     '_'.join([fail,caso,'acum{:02d}'.format(acum),fig_type]))+fmt, dpi=150, bbox_inches='tight')
        # HAGO LA LEGEND APARTE 
        # extraido de https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib
        handles, labels = ax.get_legend_handles_labels()
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=True)
        plt.close()
        def export_legend(legend, filename, expand=[-5,-5,5,5],ruta_figs=ruta_figs):
            fig = legend.figure
            fig.canvas.draw()
            bbox = legend.get_window_extent()
            bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
            bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(os.path.join(ruta_figs,filename), dpi='figure', bbox_inches=bbox)
            plt.close()
        export_legend(legend, filename='_'.join([fail,'legend',fig_type]))


    elif fig_type == 'split':
        fig, axes = plt.subplots(nrows=2,figsize=(10,10),sharex=True)
        for ax in axes.ravel():
            for o, obs in enumerate(observaciones):
                sli = slice(time_obs[o].index(time_common[2]), time_obs[o].index(time_common[-1]))
                ax.plot(time_obs[o][sli], ppmean_obs[o][sli],
                        color='black',
                        linestyle=linestyles[o],
                        linewidth=2.,
                        label=obs.name)

        for c in range(round(len([*wrfs])/len(configs))):
            estas_sims = [*wrfs][len(configs)*c:len(configs)*(c+1)]
            sim_name = estas_sims[0][1]
            lower = minimo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
            upper = maximo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
            all_array = np.zeros( (len(estas_sims), len(time_common)) )
            for s, sim in enumerate(estas_sims):
                all_array[s,:] = ppmean_wrf[sim]
            media = np.mean(all_array, axis=0)
            if wrfs[estas_sims[0]].size == 'HR':
                ax = axes.ravel()[0]
            else:
                ax = axes.ravel()[1]
            ax.fill_between(time_wrf[0][2:], lower[2:] , upper[2:],
                    color=wrfs[estas_sims[0]].color,
                    linestyle=wrfs[estas_sims[0]].linestyle,
                    linewidth=2,
                    #label=wrfs[estas_sims[0]].label,
                    alpha=0.4,
                    )
            ax.plot(time_wrf[0][2:], media[2:],
                    color=wrfs[estas_sims[0]].color,
                    linestyle=wrfs[estas_sims[0]].linestyle,
                    linewidth=2,
                    label=wrfs[estas_sims[0]].label,
                    )
        
        yticks = np.arange(0,14,2)
        for subplot,ax in enumerate(axes.ravel()):
            ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
            ax.set_yticks(list(yticks))
            ax.set_yticklabels(
                    ['{}'.format(tick) for tick in ax.get_yticks()])
            ax.set_ylabel('mm')
            ax.set_xticklabels([time[:-6] for time in time_common[2:]], fontsize=14)
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
            # HAGO LA LEGEND APARTE
            ax.grid()
            ax.set_title(casito, fontsize=14) if subplot == 0 else None
        #%% SAVE FIG
        #plt.show()
        for fmt in ['.jpg']:
            fig.savefig(os.path.join(ruta_figs,
                                     '_'.join([fail,caso,'acum{:02d}'.format(acum),fig_type]))+fmt, dpi=150, bbox_inches='tight')
    elif fig_type == 'todasSeparadas':
        cantidad_simulaciones = round(len([*wrfs])/len(configs))
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10*3,8*3),sharex=True, sharey=True)
        for ax in axes.ravel():
            for o, obs in enumerate(observaciones):
                sli = slice(time_obs[o].index(time_common[2]), time_obs[o].index(time_common[-1]))
                ax.plot(time_obs[o][sli], ppmean_obs[o][sli],
                        color='black',
                        linestyle=linestyles[o],
                        linewidth=2.,
                        label=obs.name)
        for c in range(cantidad_simulaciones):
            estas_sims = [*wrfs][len(configs)*c:len(configs)*(c+1)]
            sim_name = estas_sims[0][1]
            lower = minimo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
            upper = maximo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
            all_array = np.zeros( (len(estas_sims), len(time_common)) )
            for s, sim in enumerate(estas_sims):
                all_array[s,:] = ppmean_wrf[sim]
            media = np.mean(all_array, axis=0)
            ax = axes.ravel()[c]
            ax.fill_between(time_wrf[0][2:], lower[2:], upper[2:],
                    color=wrfs[estas_sims[0]].color,
                    linestyle=wrfs[estas_sims[0]].linestyle,
                    linewidth=2,
                    #label=wrfs[estas_sims[0]].label,
                    alpha=0.4,
                    )
            ax.plot(time_wrf[0][2:], media[2:],
                    color=wrfs[estas_sims[0]].color,
                    linestyle=wrfs[estas_sims[0]].linestyle,
                    linewidth=2,
                    label=wrfs[estas_sims[0]].label,
                    )
        yticks = np.arange(0,14,2)
        for subplot,ax in enumerate(axes.ravel()):
            ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
            ax.set_yticks(list(yticks))
            ax.set_yticklabels(
                    ['{}'.format(tick) for tick in ax.get_yticks()])
            ax.set_ylabel('mm')
            ax.set_xticklabels([time[:-6] for time in time_common[2:]], fontsize=14)
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
            # HAGO LA LEGEND APARTE
            ax.grid()
            ax.set_title(casito, fontsize=14) if subplot == 0 else None
        #%% SAVE FIG
        #plt.show()
        for fmt in ['.jpg']:
            fig.savefig(os.path.join(ruta_figs,
                '_'.join([fail,caso,'acum{:02d}'.format(acum),fig_type]))+fmt, dpi=150, bbox_inches='tight')
    elif fig_type == 'splitNestCP':
        fig, axes = plt.subplots(nrows=2,figsize=(10,10),sharex=True)
        for ax in axes.ravel():
            for o, obs in enumerate(observaciones):
                sli = slice(time_obs[o].index(time_common[2]), time_obs[o].index(time_common[-1]))
                ax.plot(time_obs[o][sli], ppmean_obs[o][sli],
                        color='black',
                        linestyle=linestyles[o],
                        linewidth=2.,
                        label=obs.name)
        for c in range(round(len([*wrfs])/len(configs))):
            estas_sims = [*wrfs][len(configs)*c:len(configs)*(c+1)]
            if wrfs[estas_sims[0]].size == 'HR':
                sim_name = estas_sims[0][1]
                lower = minimo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
                upper = maximo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
                all_array = np.zeros( (len(estas_sims), len(time_common)) )
                for s, sim in enumerate(estas_sims):
                    all_array[s,:] = ppmean_wrf[sim]
                media = np.mean(all_array, axis=0)
                if not wrfs[estas_sims[0]].label.endswith('noNest'):
                    ax = axes.ravel()[0]
                elif wrfs[estas_sims[0]].label.endswith('noNest'):
                    ax = axes.ravel()[1]
                ax.fill_between(time_wrf[0][2:], lower[2:] , upper[2:],
                        color=wrfs[estas_sims[0]].color,
                        linestyle=wrfs[estas_sims[0]].linestyle,
                        linewidth=2,
                        #label=wrfs[estas_sims[0]].label,
                        alpha=0.4,
                        )
                ax.plot(time_wrf[0][2:], media[2:],
                        color=wrfs[estas_sims[0]].color,
                        linestyle=wrfs[estas_sims[0]].linestyle,
                        linewidth=2,
                        label=wrfs[estas_sims[0]].label,
                        )
                print(wrfs[estas_sims[0]].color, wrfs[estas_sims[0]].label)
            yticks = np.arange(0,14,2)
            for subplot,ax in enumerate(axes.ravel()):
                ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
                ax.set_yticks(list(yticks))
                ax.set_yticklabels(
                        ['{}'.format(tick) for tick in ax.get_yticks()])
                ax.set_ylabel('mm')
                ax.set_xticklabels([time[:-6] for time in time_common[2:]], fontsize=14)
                for tick in ax.get_xticklabels():
                    tick.set_rotation(90)
                # HAGO LA LEGEND APARTE
                ax.grid()
                ax.set_title(casito, fontsize=14) if subplot == 0 else None
        #%% SAVE FIG
        for fmt in ['.jpg']:
            fig.savefig(os.path.join(ruta_figs,
                '_'.join([fail,caso,'acum{:02d}'.format(acum),fig_type]))+fmt, dpi=150, bbox_inches='tight')
        #plt.show()

    elif fig_type == '3WaySplit':
        fig, axes = plt.subplots(nrows=3, figsize=(10,10), sharex=True)
        for ax in axes.ravel():
            for o, obs in enumerate(observaciones):
                sli = slice(time_obs_dt[o].index(time_common_dt[2]), time_obs[o].index(time_common[-1]))
                ax.plot(time_obs_dt[o][sli], ppmean_obs[o][sli],
                        color='black',
                        linestyle=linestyles[o],
                        linewidth=2,
                        label=obs.name)
        for c in range(round(len([*wrfs])/len(configs))):
            estas_sims = [*wrfs][len(configs)*c:len(configs)*(c+1)]
            if wrfs[estas_sims[0]].size == 'CP':
                sim_name = estas_sims[0][1]
                lower = minimo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
                upper = maximo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
                all_array = np.zeros( (len(estas_sims), len(time_common)) )
                for s, sim in enumerate(estas_sims):
                    all_array[s,:] = ppmean_wrf[sim]
                media = np.mean(all_array, axis=0)

                if wrfs[estas_sims[0]].label.startswith('CP1'):
                    ax = axes.ravel()[0]
                elif not wrfs[estas_sims[0]].label.endswith(')'):
                    ax = axes.ravel()[2]
                else:
                    ax = axes.ravel()[1]
            
                ax.fill_between(time_wrf_dt[0][2:], lower[2:], upper[2:],
                        color=wrfs[estas_sims[0]].color,
                        linestyle=wrfs[estas_sims[0]].linestyle,
                        linewidth=2,
                        alpha=0.4,
                        )
                ax.plot(time_wrf_dt[0][2:], media[2:],
                        color=wrfs[estas_sims[0]].color,
                        linestyle=wrfs[estas_sims[0]].linestyle,
                        linewidth=2,
                        label=wrfs[estas_sims[0]].label,
                        )
            yticks = np.arange(0,14,2)
            for subplot,ax in enumerate(axes.ravel()):
                ax.set_ylim(np.min(yticks)-1, np.max(yticks)+1)
                ax.set_yticks(list(yticks))
                ax.set_yticklabels(
                        ['{}'.format(tick) for tick in ax.get_yticks()])
                ax.set_ylabel('mm')
                ax.set_xlim(time_common_dt[2]-datetime.timedelta(hours=6), time_common_dt[-1]+datetime.timedelta(hours=0))
                ax.xaxis.set_major_locator(mdates.DayLocator())
                ax.xaxis.set_minor_locator(mdates.HourLocator((3,6,9,12,15,18,21)))
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%d'))
                ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
                ax.xaxis.set_tick_params(which='major', pad=16) 
                ax.grid()
                ax.grid(which='minor', alpha=0.4)
                ax.set_title(name_casito[casito], fontsize=20) if subplot == 0 else None
        fig.subplots_adjust(hspace=0.1)
        #%% SAVE FIG
        for fmt in ['.png', '.pdf']:
            fig.savefig(os.path.join(ruta_figs,
                '_'.join([fail,caso,'acum{:02d}'.format(acum),fig_type]))+fmt, dpi=150, bbox_inches='tight')
        #plt.show()
        


        ### GUARDO LAS LEGENDS APARTE ###
        fig, axes = plt.subplots(figsize=(10,10))
        for o, obs in enumerate(observaciones):
            sli = slice(time_obs_dt[o].index(time_common_dt[2]), time_obs[o].index(time_common[-1]))
            ax.plot(time_obs_dt[o][sli], ppmean_obs[o][sli],
                    color='black',
                    linestyle=linestyles[o],
                    linewidth=2,
                    label=obs.name)

        for c in range(round(len([*wrfs])/len(configs))):
            estas_sims = [*wrfs][len(configs)*c:len(configs)*(c+1)]
            if wrfs[estas_sims[0]].size == 'CP':
                sim_name = estas_sims[0][1]
                lower = minimo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
                upper = maximo_de_varios( [ppmean_wrf[j] for j in estas_sims ] )
                all_array = np.zeros( (len(estas_sims), len(time_common)) )
                for s, sim in enumerate(estas_sims):
                    all_array[s,:] = ppmean_wrf[sim]
                media = np.mean(all_array, axis=0)
                ax.fill_between(time_wrf_dt[0][2:], lower[2:], upper[2:],
                        color=wrfs[estas_sims[0]].color,
                        linestyle=wrfs[estas_sims[0]].linestyle,
                        linewidth=2,
                        alpha=0.4,
                        )
                ax.plot(time_wrf_dt[0][2:], media[2:],
                        color=wrfs[estas_sims[0]].color,
                        linestyle=wrfs[estas_sims[0]].linestyle,
                        linewidth=2,
                        label=wrfs[estas_sims[0]].label,
                        )
        
        handles, labels = ax.get_legend_handles_labels()
        legend = plt.legend(handles, labels, ncol=round(len([*wrfs])/len(configs))+3, loc=3, framealpha=1, frameon=True, bbox_to_anchor=(1.05,1))
        def export_legend(legend, filename=os.path.join(ruta_figs,'_'.join([fail,caso,'acum{:02d}'.format(acum),fig_type, 'legend']))+'.png'):
            fig  = legend.figure
            fig.canvas.draw()
            bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(filename, dpi="figure", bbox_inches=bbox)

        export_legend(legend)

