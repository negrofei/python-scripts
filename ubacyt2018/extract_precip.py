"""

 ### UBACYT2018 ###

 Script para extraer a un archivo netcdf la precipitación acumulada en determinado periodo "acum" para corridas de WRF ARW. 

ejecutar como:
    $ python3 extrac_precip.py <wrf_netcdf_file> 

devuelve un archivo llamado precip_<acum>_<wrf_netcdf_file>.nc

"""

from netCDF4 import Dataset, date2num
import Docpy
import numpy as np
import time as time
from datetime import datetime, timedelta
import os
import sys

# INPUT
ruta = os.path.abspath(sys.argv[1])     
acum = '{:02d}'.format(int(input('Elegí acumulacion (hs):\t')))

wrfori = Docpy.WRF(ruta)
new_wrf_name = '_'.join(['precip',acum,wrfori.filename])

str_tiempos, delta_t_hs = wrfori.timedata()
dataset = Dataset(os.path.join(os.path.split(ruta)[0],new_wrf_name), 'w', format=wrfori.data_netcdf.file_format)
# Traigo las variables del original

lats, lons = wrfori.get_latlon()
# Creo las dimensiones del archivo en base al archivo original 
dataset.createDimension('lat', len(lats))
dataset.createDimension('lon', len(lons))
dataset.createDimension('time', None)

# Creo las variables con las dimensiones 
# Tiempos
times = dataset.createVariable('time', np.float64, ('time',))
times.units = 'hours since 1-1-1 00:00:00'
times.calendar = 'standard'

# Latitudes
latitude = dataset.createVariable('lat', np.float64, ('lat',))
latitude[:] = lats
latitude.units = 'degree_north'

# Longitudes
longitude = dataset.createVariable('lon', np.float64, ('lon',))
longitude[:] = lons
longitude.units = 'degree_east'

# Información del archivo
dataset.description = ' '.join(['Precipitacion extraida desde',
                                os.path.abspath(ruta),
                                'a partir del script:',
                                os.path.abspath(sys.argv[0])])
dataset.history = '\n'.join(['Archivo creado el: ' + time.ctime( time.time() ), wrfori.data_netcdf.history])
dataset.source = ''

# Create the actual 4-d variable
precip = dataset.createVariable('precip', np.float64, ('time','lat','lon'))
precip.units = 'mm'
precip.long_name = 'Precipitación acumulada HACIA ADELANTE en {} hs '.format(acum)

precip[:] = Docpy.precip.calc_precip_acum(wrfori, acum=int(acum))
times[:] = [date2num(datetime.strptime(Docpy.functions.acum_time(wrfori,acum=int(acum))[j],'%Y-%m-%d %H:%M:%S'), units = times.units, calendar = times.calendar) for j in range(len(Docpy.functions.acum_time(wrfori, acum=int(acum))))]

dataset.close()
print('Archivo',new_wrf_name,'creado!')

