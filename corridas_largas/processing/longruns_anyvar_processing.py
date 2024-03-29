"""
This script takes the directory of a long simulation with tesla-larga as argument
and returns the wrfouts processed divided in files per field

for execution

$ python3 longruns_anyvar_processing.py <dir_of_chunks>

returns

<var_name>_every**_alltimes_wrfout_<dir_of_chunks>.nc

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
import concurrent.futures
import itertools
from pathlib import Path
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Defino la función que procesa un archivo

def process_wrfouts( wrfouts_path, var, levels=[1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.,250.,200.,150.,100.]):
    #print('PROCESS:', os.getpid())
    
    chunk = Path(wrfouts_path[0]).parent
    print('Chunk:',chunk.name, var)
    #print('Leo los datasets')
    wrflist = [Dataset(wrfout) for wrfout in wrfouts_path]
    wrfori = wrflist[0]     # agarro un wrf para sacar metadata
    wrffin = wrflist[-1]
    #print('Obtengo tiempos')
    tiempos =  wrf.getvar(wrflist, 'XTIME', timeidx=wrf.ALL_TIMES, method='cat') 
    time_list = [str(num2date(tiempos[i], units=tiempos.units)) for i in range(len(tiempos[:]))]
    delta_t_hs = num2date(tiempos[1]-tiempos[0], units=tiempos.units).hour
    every = '03'
    #print('Outout every ',every,'hours')
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
    
    
    new_wrf_name_beginning = filename_generator(var, levels)

    new_wrf_name = '_'.join([new_wrf_name_beginning,'every', every, wrf_chunk_name])
    #print('Nombre del archivo:',new_wrf_name)

    dataset = Dataset( 
            chunk.joinpath(new_wrf_name),
            'w',
            format=wrfori.file_format
            )
    #print('Creating lonlats...')
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

    #print('Creating time dimensions...')
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
    #print('Creating file info...')
    dataset.description = ' '.join(
            [
                'Geopotencial extraido desde',
                str(chunk),
                'a partir del script:',
                os.path.abspath(sys.argv[0]),
                ],
            )

    dataset.history = '\n'.join(['Archivo creado el: ' + time.ctime( time.time() ), ])
    dataset.source = ''

    ### Extraigo la variable
    # Primero chequeo si es de 2 o 3 dimensions
    shape = wrf.getvar(wrflist, var, timeidx=0, method='cat').shape
    if np.size(shape) == 3:
        array = np.zeros( (len(tiempos), z,y,x) )
        dims = tuple(['time','lev','lat','lon'])
        for t in range(len(time_list[::inter])):
            array[t,:] = wrf.interplevel(
                                         field3d=wrf.getvar( wrflist, var, timeidx=t, method='cat'),
                                         vert=p,
                                         desiredlev=levels)
    elif np.size(shape) == 2:
        array = np.zeros( (len(tiempos), 1,y,x) )
        dims = tuple(['time','lev_2','lat','lon'])
        for t in range(len(time_list[::inter])):
            array[t,:] = wrf.getvar( wrflist, var, timeidx=t, method='cat' )

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
    return new_wrf_name


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def check_if_var_exists(variable):
    assert corrida_larga_name
    
    dir_to_mokdata = sorted(                                # Esto me deja con el primer directorio
            os.listdir(                                     # para levantar un archivo cualquiera
                os.path.abspath(corrida_larga_name)
                )
            )[0]

    file_to_mokdata = list(                                     # Esto me agarra los archivos en el 
            filter(                                             # dir_to_mok y se queda solo con los
                lambda x: x.startswith('wrfout_d0'), sorted(    # que empiezan con wrfout_d0
                    os.listdir(                                 # luego me quedo con el primer item
                        os.path.join(corrida_larga_name, dir_to_mokdata)    
                        )
                    )
                )
            )[0]

    # Chequeo que la variable sea válida
    mokdata = Dataset(
            os.path.join(
                './',
                corrida_larga_name,
                dir_to_mokdata,
                file_to_mokdata))
            
    assert (variable in [wrf.routines._FUNC_MAP]) or (variable in [*mokdata.variables]) or (wrf.routines._undo_alias(variable) in wrf.routines._VALID_KARGS), "{} should be of the available wrf-python diagnostics or a WRFOUT valid variable".format(variable)
    return mokdata


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def filename_generator(variable, levels=None):
 
    mokdata = check_if_var_exists(variable)
    shape = wrf.getvar( mokdata, variable).shape
    if np.size(shape) == 2:
        lev_name = 'single_lvl'
    elif np.size(shape) == 3:
        if levels:
            # El nombre del archivo tiene que contener el nivel o cantidad de niveles verticales.
            if len(levels)==1:
                lev_name = str(int(levels[0]))+'hpa'
            elif len(levels)>1:
                lev_name = str(len(levels))+'lvls'

    new_wrf_name_beginning = '_'.join([variable,lev_name])
    return new_wrf_name_beginning
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Funcion para ensamblar los archivos con CDO

def concatenate( list_of_files, dirout, nameout ):
    list_of_files = [str(a) for a in list_of_files]
    c = cdo.Cdo()
    c.cat(input=list_of_files, output=os.path.join(dirout, nameout))
    print(nameout)
    return nameout 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Funcion para dar los wrfouts_paths segun el chunk

def wrfpaths(chunk, corrida_larga_name):

    wrfouts_path = sorted([corrida_larga_name.joinpath(chunk, file) for file in corrida_larga_name.joinpath(chunk).iterdir() if file.name.startswith('wrfout')])
    return wrfouts_path

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    tini = time.perf_counter()
    ### Directorios

    corrida_larga_name = Path(sys.argv[1]).resolve()                        # Carpeta donde se encuentran los files
    #variable = input('Choose the variable to extract:\t')   # Variable para procesar

    # Voy a automatizar el proceso de seleccion de variables con la siguiente lista
    # (RAINS Y MCAPE DEBEN HACERSE POR SEPARADO EN OTRO SCRIPT):
    #
    #   *   avo         |   Absolute Vorticity
    #   *   eth         |   Equivalent Potential Temperature
    #   *   mdbz        |   Maximum Reflectivity
    #   *   geopt       |   Geopotential for the Mass Grid
    #   *   omg         |   Omega
    #   *   pressure    |   Full Model Pressure
    #   *   rh          |   Relative Humidity
    #   *   tc          |   Temperature in Celcius
    #   *   td          |   Dew Point Temperature
    #   *   ua          |   U-component of Wind on Mass Points
    #   *   va          |   V-component of Wind on Mass Points
    #   *   wa          |   W-component of Wind on Mass Points
    #   *   QVAPOR      |   Water vapor mixing ratio
    #   *   OLR         |   TOA OUTGOING LONGWAVE

    variables = ['T2','avo','eth']#,'mdbz', 'geopt','omg','pressure','rh','tc','td','ua','va','wa','QVAPOR','OLR']

    ### Voy a crear hacer threading de forma de que haga cada variable y cada chunk en un proceso 
    #   aparte y que no espere que termine para lanzar el otro proceso

    # Tengo que iterar entre los chunks y entre las variables. 
   
    # Agarro los chunks
    chunks = sorted(corrida_larga_name.glob('20*'))
    
    # Obtengo la lista de parametros para pasarle al threading
    parametros = list(itertools.product(chunks, variables))
    
    # Como tengo que mapear process_wrfouts, armo una funcion auxiliar para parsear los parametros
    def parser_wrfprocess(args):
        return process_wrfouts( wrfpaths(args[0], corrida_larga_name), args[1] )
    
    ### Una vez termina de armar el <variable>_every*_<chunk>.nc lo ensamblo en un solo archivo
    #   Como puse los threads, lo que va a pasar es que hasta que no se termine de ejecutar
    #   eso no puede pasar a lo siguiente. Esto no es lo más optimo. 
    #   Lo más óptimo sería que la siguiente parte se ejecute al mismo tiempo y que chequee que 
    #   estén todos los archivos. Armo una funcion que haga eso. 

    def cdo_ensemble(variable, chunks):
        # Chequeo que esten todos los chunks terminados
        list_of_files = []
        while True:
            for chunk in chunks:
                file = next(corrida_larga_name.joinpath(chunk).glob(variable+'*'))
                list_of_files = list_of_files + [file]

            if len(list_of_files) == len(chunks):
                printer('Starting to ensemble', variable)
                break
            else:
                print('Waiting for', variable,'\t','...','{}/{} files done'.format(len(list_fo_files),len(chunks)))
                time.sleep(2)

        # Me quedo con la raiz del nombre de los archivos, ya que contiene la info de cuántos niveles
        # hay
        printer('ENMSABLE OF',variable)

        root = list_of_files[0].name.split('_wrfout')[0]

        dirout = str(corrida_larga_name.resolve())
        nameout = '_'.join(
                [
                    root,
                    'alltimes',
                    'wrfout',
                    corrida_larga_name.name+'.nc'
                    ]
                )
        
        concatenate(list_of_files, dirout, nameout)
        return nameout
    
    def parser_ensemble(args):
        return cdo_ensemble(args, chunks)


    # Arranco el threading
    with concurrent.futures.ThreadPoolExecutor(20) as executor:
        results_wrfprocess = executor.map(parser_wrfprocess, parametros)
    
    results_ensamble = map(parser_ensemble, variables)

    # Performance    
    tfin = time.perf_counter()

    printer('Tardo','{:.1f}'.format(tfin-tini))

