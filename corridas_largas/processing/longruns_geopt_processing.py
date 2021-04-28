"""
This script takes the directory of a long simulation with tesla-larga as argument
and returns the wrfouts processed divided in files per field

for execution

$ python3 longruns_geopt_processing.py <dir_of_chunks>

returns

geopt_every**_alltimes_wrfout_<dir_of_chunks>.nc

"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Paquetes

from Docpy.functions import printer
import numpy as np
import os 
import sys
from netCDF4 import Dataset, num2date, date2num
import wrf
import time
from datetime import datetime, timedelta
import cdo
import glob

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Defino la función que procesa un archivo

def process_wrfouts( wrfouts_path, levels=[1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.,250.,200.,150.,100.]):
    print('PROCESS:', os.getpid())
    
    chunk = os.path.split(wrfouts_path[0])[0]
    print('Chunk:',chunk)
    print('Leo los datasets')
    wrflist = [Dataset(wrfout) for wrfout in wrfouts_path]
    wrfori = wrflist[0]     # agarro un wrf para sacar metadata
    wrffin = wrflist[-1]
    print('Obtengo tiempos')
    tiempos =  wrf.getvar(wrflist, 'XTIME', timeidx=wrf.ALL_TIMES, method='cat') 
    time_list = [str(num2date(tiempos[i], units=tiempos.units)) for i in range(len(tiempos[:]))]
    delta_t_hs = num2date(tiempos[1]-tiempos[0], units=tiempos.units).hour
    every = '03'
    print('Outout every ',every,'hours')
    inter = int(int(every)/delta_t_hs)
    
    ### Extraigo partes del código de script de extract_precip. 
    wrf_chunk_name = '_'.join(
            [
                os.path.split(wrfori.filepath())[-1][:-9],  
                # me quedo con el string 'wrfout_d0*_AAAA-MM-DD
                os.path.split(wrffin.filepath())[-1].split('_')[2]+'.nc',  
                # me quedo con la fecha de fin del chunk AAAA-MM-DD
                ]
            )
    # El nombre del archivo tiene que contener el nivel o cantidad de niveles verticales.

    if len(levels)==1:
        lev_name = str(int(levels[0]))+'hpa'
    elif len(levels)>1:
        lev_name = str(len(levels))+'lvls'

    new_wrf_name = '_'.join(['geopt_every',every, lev_name, wrf_chunk_name])
    print('Nombre del archivo:',new_wrf_name)

    dataset = Dataset( 
            os.path.join(chunk, new_wrf_name),
            'w',
            format=wrfori.file_format
            )
    print('Creating lonlats...')
    lats = wrf.getvar(wrfori, 'lat')
    lons = wrf.getvar(wrfori, 'lon')

    # Creo las dimensiones del archivo en base al archivo original 
    dataset.createDimension('lev', len(levels))
    dataset.createDimension('lev_2', 1)
    dataset.createDimension('lat', lats.shape[0])
    dataset.createDimension('lon', lons.shape[1])
    dataset.createDimension('time', None)

    # Voy a interpolar a niveles de presión
    p = wrf.getvar( wrflist, 'pressure' )
    z = len(levels)
    y, x = p.shape[1:]

    print('Creating time dimensions...')
    # Creo las variables con las dimensiones 
    # Tiempos
    times = dataset.createVariable('time', np.float64, ('time',))
    times.units = 'hours since 1-1-1 00:00:00'
    times.calendar = 'standard'

    # Niveles
    levs = dataset.createVariable('lev', np.float64, ('lev',))
    levs[:] = levels
    levs.units = 'pressure (hPa)'

    #Niveles2
    levs2 = dataset.createVariable('lev_2', np.float64, ('lev_2',))
    levs2[:] = 1000.0,
    levs2.units = 'hPa'

    # Latitudes
    latitude = dataset.createVariable('lat', np.float64, ('lat',))
    latitude[:] = lats[:,0]
    latitude.units = 'degree_north'

    # Longitudes
    longitude = dataset.createVariable('lon', np.float64, ('lon',))
    longitude[:] = lons[0,:]
    longitude.units = 'degree_east'

    # Información del archivo
    print('Creating file info...')
    dataset.description = ' '.join(
            [
                'Geopotencial extraido desde',
                os.path.abspath(chunk),
                'a partir del script:',
                os.path.abspath(sys.argv[0]),
                ],
            )

    dataset.history = '\n'.join(['Archivo creado el: ' + time.ctime( time.time() ), ])
    dataset.source = ''

    # Extraigo la geopt 
    varis = ['geopt']
    for var in varis:
        
        array = np.zeros( (len(tiempos), z,y,x) )
        for t in range(len(time_list[::inter])):
            array[t,:] = wrf.interplevel(
                                         field3d=wrf.getvar( wrflist, var, timeidx=t, method='cat'),
                                         vert=p,
                                         desiredlev=levels)

        dims = tuple(['time','lev','lat','lon'])
        variable = dataset.createVariable(var.lower(), np.float64, dims)
        # unidades y description
        campo = wrf.getvar(wrflist, var, timeidx=t, method='cat')
        variable.units = campo.units
        variable.description = campo.description
        variable[:] = array

    times[:] = [
            date2num(
                datetime.strptime(
                    time_list[::inter][j],
                    '%Y-%m-%d %H:%M:%S'
                    ),
                units = times.units,
                calendar = times.calendar
                )
            for j in range(len(time_list[::inter]))
            ]

    dataset.close()
    print('Archivo',new_wrf_name,'creado!')
    return None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Funcion para ensamblar los archivos con CDO

def concatenate( list_of_files, dirout, nameout ):
    c = cdo.Cdo()
    c.cat(input=list_of_files, output=os.path.join(dirout,nameout) )
    return None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':
    ### Directorios

    corrida_larga_name = sys.argv[1]
    chunks = glob.glob(os.path.join(os.path.abspath(corrida_larga_name),'20*'))
    chunks.sort()
    for chunk in chunks:
        wrfouts_path = [os.path.join(corrida_larga_name, chunk,file) for file in os.listdir(os.path.join(corrida_larga_name, chunk)) if file.startswith('wrfout')]
        wrfouts_path.sort()

        process_wrfouts( wrfouts_path, levels=[850.] )
    
    ### una vez termina de armar el geopt_every*_<chunk>.nc lo ensamblo en un solo archivo
    list_of_files = []
    for chunk in chunks:
        file = glob.glob(os.path.join(corrida_larga_name,chunk,'geopt*'))
        list_of_files = list_of_files + file
        print(file)
    
    # Me quedo con la raiz del nombre de los archivos, ya que contiene la info de cuántos niveles
    # hay
    
    root = os.path.split(file[0])[-1].split('_wrfout')[0]

    dirout = os.path.abspath(corrida_larga_name)
    nameout = '_'.join(
            [
                root,
                'alltimes',
                'wrfout',
                corrida_larga_name+'.nc'
                ]
            )
    
    concatenate(list_of_files, dirout, nameout)

