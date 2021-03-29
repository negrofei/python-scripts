"""
This script takes the directory of a long simulation with tesla-larga as argument
and returns the wrfouts processed divided in files per field
"""

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
# defino la función que procesa un archivo

def process_wrfouts( wrfouts_path ):
    print('PROCESS:', os.getpid())
    
    chunk = os.path.split(wrfouts_path[0])[0]
    print('Chunk:',chunk)
    print('Leo los datasets')
    wrflist = [Dataset(wrfout) for wrfout in wrfouts_path]
    wrfori = wrflist[0] # agarro un wrf para sacar metadata
    wrffin = wrflist[-1]
    print('Obtengo tiempos')
    tiempos =  wrf.getvar(wrflist, 'XTIME', timeidx=wrf.ALL_TIMES, method='cat') 
    time_list = [str(num2date(tiempos[i], units=tiempos.units)) for i in range(len(tiempos[:]))]
    delta_t_hs = num2date(tiempos[1]-tiempos[0], units=tiempos.units).hour
    acum = '{:02d}'.format(3)
    print('Acumulo cada ',acum,'horas') 
    inter = int(int(acum)/delta_t_hs)
#    print('Obtengo la RAIN')    
#    rain = wrf.getvar(wrflist, 'RAINC', timeidx=wrf.ALL_TIMES, method='cat') + wrf.getvar(wrflist, 'RAINSH', timeidx=wrf.ALL_TIMES, method='cat') + wrf.getvar(wrflist, 'RAINNC', timeidx=wrf.ALL_TIMES, method='cat')
    
    ### Extraigo partes del código de script de extract_precip. 
    wrfori_name = '_'.join([
        os.path.split(wrfori.filepath())[-1][:-9],
        os.path.split(wrffin.filepath())[-1].split('_')[2]+'.nc'
        ])
    
    new_wrf_name = '_'.join(['precip_acum',acum, wrfori_name])
    print('Nombre del archivo:',new_wrf_name)

    dataset = Dataset( 
            os.path.join(chunk,new_wrf_name),
            'w',
            format=wrfori.file_format
            )
    print('Creating lonlats...')
    lats = wrf.getvar(wrflist[0], 'lat')
    lons = wrf.getvar(wrflist[0], 'lon')

    # Creo las dimensiones del archivo en base al archivo original 
    dataset.createDimension('lev_2',1)
    dataset.createDimension('lat', lats.shape[0])
    dataset.createDimension('lon', lons.shape[1])
    dataset.createDimension('time', None)

    print('Creating time dimensions...')
    # Creo las variables con las dimensiones 
    # Tiempos
    times = dataset.createVariable('time', np.float64, ('time',))
    times.units = 'hours since 1-1-1 00:00:00'
    times.calendar = 'standard'

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
    dataset.description = ' '.join(['Precipitacion extraida desde',
                                    os.path.abspath(chunk),
                                    'a partir del script:',
                                    os.path.abspath(sys.argv[0])])

    dataset.history = '\n'.join(['Archivo creado el: ' + time.ctime( time.time() ), ])
    dataset.source = ''

#    # Creo la variable de precipitacion
#    print('Creating precipitation acum',acum)
#    precip = dataset.createVariable('precip', np.float64, ('time','lat','lon'))
#    precip.units = 'mm'
#    precip.long_name = 'Precipitación acumulada HACIA ATRAS en {} hs '.format(acum)
#    
#    inter = int(int(acum)/delta_t_hs)
#    precip_acum = np.ma.empty([len(time_list[::inter]),rain.shape[1],rain.shape[2]])
#    for r in range(1,precip_acum.shape[0]):
#        precip_acum[r,:,:] = rain[(r)*inter,:]-rain[r-1*inter,:]

    # No voy a acumular la precipitacion acá, sino voy a hacerlo compatible con Docpy
    rainvars = ['RAINNC','RAINC','RAINSH']
    for var in rainvars:
        shape = wrf.getvar(wrflist, var, timeidx=0, method='cat').shape
        array = np.zeros( (len(tiempos), 1,)+shape )
        for t in range(len(time_list[::inter])):
            array[t,:] = wrf.getvar( wrflist, var, timeidx=t, method='cat')

        dims = tuple(['time','lev_2','lat','lon'])
        variable = dataset.createVariable(var.lower(), np.float64, dims)
        # unidades y description
        campo = wrf.getvar(wrflist, var, timeidx=t, method='cat')
        variable.units = campo.units
        variable.description = campo.description
        variable[:] = array

    times[:] = [date2num(datetime.strptime(time_list[::inter][j],'%Y-%m-%d %H:%M:%S'), units = times.units, calendar = times.calendar) for j in range(len(time_list[::inter]))]

    dataset.close()
    print('Archivo',new_wrf_name,'creado!')
    return None

# funcion para ensamblar los archivos con CDO
def concatenate( list_of_files, dirout, nameout ):
    c = cdo.Cdo()
    c.cat(input=list_of_files, output=os.path.join(dirout,nameout) )
    return None

if __name__ == '__main__':
    ### Directorios

    corrida_larga_name = sys.argv[1]
    chunks = os.listdir(os.path.abspath(corrida_larga_name))
    chunks.sort()
    for chunk in chunks:
        wrfouts_path = [os.path.join(corrida_larga_name, chunk,file) for file in os.listdir(os.path.join(corrida_larga_name, chunk)) if file.startswith('wrfout')]
        wrfouts_path.sort()

        process_wrfouts( wrfouts_path )
    
    ### una vez termina de armar el precip_acum03.... lo ensamblo en un solo archivo
    list_of_files = []
    for chunk in chunks:
        file = glob.glob(os.path.join(corrida_larga_name,chunk,'precip*'))
        list_of_files = list_of_files + file
        print(file)
    dirout = os.path.abspath(corrida_larga_name)
    nameout = '_'.join([ \
                        'precip',\
                        'acum03',\
                        'alltimes',\
                        'wrfout',\
                        corrida_larga_name+'.nc'])
    
    concatenate(list_of_files, dirout, nameout)

