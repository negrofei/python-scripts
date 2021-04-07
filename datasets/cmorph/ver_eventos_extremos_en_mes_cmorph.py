"""
Script para ver eventos de precipitación diara en un mes determinado en una caja determinada.


"""
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Paquetes 
import os
import sys
import glob
import re
import Docpy
from Docpy import Dato

import numpy as np

import matplotlib.pyplot as plt

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Leo los archivos mensuales que empiezan con CMORPH y tienen un solo mes

archivos = glob.glob('./CMORPH*')
expr = re.compile(r'(_\d{6}_)')
archivos = [archivo for archivo in archivos if expr.search(archivo)]
archivos.sort()

datas = dict(zip(archivos, [Dato(archivo) for archivo in archivos]))

### Promedio la precipitacion en la caja box4 de toda la vida

lat, lon = datas[archivos[0]].get_latlon()

box4 = Docpy.functions.utils['box']
print('Precipitation average taken in\n', box4)

lonmin_idx = np.where(np.abs(lon-box4[0])==np.min(np.abs(lon-box4[0])))
lonmax_idx = np.where(np.abs(lon-box4[1])==np.min(np.abs(lon-box4[1])))
latmin_idx = np.where(np.abs(lat-box4[2])==np.min(np.abs(lat-box4[2])))
latmax_idx = np.where(np.abs(lat-box4[3])==np.min(np.abs(lat-box4[3])))

### Calculo la precip diaria en la caja para cada mes

meses = [re.search(r'\d{6}', archivo).group() for archivo in archivos]
times = [Docpy.functions.acum_time(datas[archivos[i]], 24) for i in range(len(archivos))]
mean_pp_box = [Docpy.precip.ppmean_box(datas[archivos[i]], 24, box4) for i in range(len(archivos))]


casos = ['2005-03-12', '2015-11-10','2016-10-24']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Grafico

fig, axes = plt.subplots(figsize=(10,10), ncols=3, nrows=3, sharey=True)
newaxes = axes.reshape(-1)
for m in range(len(meses)):
    newaxes[m].plot(times[m],mean_pp_box[m], linewidth=2)
    newaxes[m].set_title(meses[m])
    newaxes[m].tick_params(axis='x', bottom=False, labelbottom=False)
    if meses[m] == '200503':
        newaxes[m].scatter(times[m][11], mean_pp_box[m][11],color='tab:red')
    elif meses[m] == '201511':
        newaxes[m].scatter(times[m][9], mean_pp_box[m][9],color='tab:red')
    elif meses[m] == '201610':
        newaxes[m].scatter(times[m][23], mean_pp_box[m][23],color='tab:red')
plt.show()

### Savefig
fail = sys.argv[0].split('.py')[0]
sape = '_'.join([fail])
ruta_figs = os.path.join('./','figuras',fail)
try:
    os.makedirs(ruta_figs)
except:
    print(ruta_figs,'ya existe')

fig.savefig(os.path.join(ruta_figs,sape), dpi=150, bbox_inches='tight')
