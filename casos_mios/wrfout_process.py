"""

Script para reemplazar el ARWpost por el modulo de wrf-python

Al final no lo terminé usando porque había algun tema de incompatibilidades que implicaba mucho retoque del paquete Docpy

"""


import time
import wrf
import glob
import os
import sys
from Docpy.functions import printer
import numpy as np
import datetime
from netCDF4 import Dataset, date2num, num2date

def process_wrfout( wrfout_list,\
        varlist=['geopt','tc','ua','va','HGT','RAINNC','RAINC','RAINSH','QVAPOR'],
        levels=[1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.,250.,200.,150.,100.],
        **kwargs):


    print( 'PROCESS:', os.getpid() )
    print( 'varlist to process:',varlist )
    datalist = [Dataset(wrfout) for wrfout in wrfout_list]
    ruta_wrfouts = os.path.split(os.path.abspath(wrfout_list[0]))[0]

    
    # cada cuanto quiero que esté el wrfout
    tiempos = wrf.getvar(datalist, 'xtimes', timeidx=wrf.ALL_TIMES)
    delta_wrfout = int( tiempos[1] -\
                        tiempos[0] )//60
    if 'every' in kwargs:
        wanted_delta = kwargs['every']
    else:
        wanted_delta = 3
    inter = wanted_delta//delta_wrfout
    if inter < 1:
        printer('WARNING')
        print('Los wrfouts están en una resolucion temporal más baja que lo que vos querés')
        print('wrfout_delta=',wrfout_delta)
        print('wanted_delta=',wanted_delta)
        print('Volviendo a la original...')
        wanted_delta = delta_wrfout
        inter= wanted_delta//delta_wrfout
    timeidxs = range(0,len(tiempos), inter)    
    
    # voy a interpolar a los niveles de presión
    p = wrf.getvar(datalist, 'pressure')
    z = len(levels)
    y,x = p.shape[1:]
    
    dimsdict = {len(tiempos):'time',
                1:'lev_2',
                z:'lev',
                y:'lat',
                x:'lon'}


    # ahora armo los netcdf4
    wrfori = datalist[0]
    namelist = wrfout_list[0].split('_')
    namelist.insert(1,'py')
    namelist[3] = namelist[3][:-3]
    namelist.pop(4)
    new_wrf_name = '_'.join(namelist)+'.nc'
    print('Nombre del archivo:',new_wrf_name)
    
    dataset = Dataset(
            os.path.join(ruta_wrfouts,new_wrf_name),
            'w',
            format=wrfori.file_format
            )
    
    # lat lon
    lats = wrf.getvar(datalist[0], 'lat')
    lons = wrf.getvar(datalist[0], 'lon')

    # Creo las dimensiones del archivo en base al archivo original 
    dataset.createDimension('lat', y)
    dataset.createDimension('lon', x)
    dataset.createDimension('lev', z)
    dataset.createDimension('lev_2',1)
    dataset.createDimension('time', None)

    print('Creating time dimensions...')
    # Creo las variables con las dimensiones 
    # Tiempos
    times = dataset.createVariable('time', np.float64, ('time',))
    times.units = 'hours since 1-1-1 00:00:00'
    times.calendar = 'standard'

    # Niveles
    levs = dataset.createVariable('lev', np.float64, ('lev',) )
    levs[:] = levels
    levs.units = 'pressure (hPa)'

    #niveles2
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
    dataset.description = ' '.join(['Archivos originales:',
                                    ruta_wrfouts,
                                    'Procesado a partir del script:',
                                    os.path.abspath(sys.argv[0])])

    dataset.history = '\n'.join(['Archivo creado el: ' + time.ctime( time.time() ), ])
    dataset.source = ''

    # creo las variables del netCDF
    # me traigo las variables
    for var in varlist:
        printer(var)
        ini = time.time()
        shape = wrf.getvar(datalist, var, timeidx=0, method='cat').shape
        if shape==(y,x):
            array = np.zeros( (len(tiempos),1,)+shape)
            for t in range(len(timeidxs)):
                array[t,:] = wrf.getvar(datalist, var, timeidx=t, method='cat')
        else:
            array = np.zeros( (len(tiempos),z,y,x) )
            for t in range(len(timeidxs)):

                array[t,:] = wrf.interplevel( 
                                            field3d=wrf.getvar(datalist, var, timeidx=t, method='cat'),
                                            vert=p,
                                            desiredlev=levels)
        dims = tuple([dimsdict[a] for a in array.shape])
        variable = dataset.createVariable(var.lower(), np.float64, dims)
        # unidades y descripcion:
        campo = wrf.getvar(datalist, var, timeidx=t, method='cat')
        variable.units = campo.units
        variable.description = campo.description
        variable[:] = array
        print( str( datetime.timedelta( seconds=time.time()-ini ) ) )
    
    # finalmente lleno la variable de tiempos
    times[:] = [date2num(num2date(tiempos[t], units=tiempos.units), units=times.units, calendar=times.calendar) for t in range(len(tiempos))]
    dataset.close()
    print('Archivo',new_wrf_name,'creado!')
    return None

if __name__ == '__main__':
    pyScript = sys.argv[0]
    printer('ACORDATE QUE EL ARGUMENTO VA ENTRE ""')
    wrfout_input = sys.argv[1]
    wrfout_list = sorted( glob.glob( wrfout_input ) )
    print('Te voy a procesesar estos archivos:')
    [print(coso) for coso in wrfout_list]
    process_wrfout( wrfout_list )
   

