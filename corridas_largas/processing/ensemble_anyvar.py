"""
This script makes the ensemble of physics of a field.

for execution:
    $ python ensemble_anyvar.py 

"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Paquetes

import cdo
import glob
import os
import pathlib
import re
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Directorios
root = pathlib.Path('.')
dir_ensemble = root.joinpath('ensemble')

# Me quedo con los directorios que contienen las físicas (arrancan con ^[Cc]\d_) y los rains

patron = re.compile(r'^[Cc]\d_')
dir_physics = [x for x in root.iterdir() if x.is_dir() and re.match(patron, str(x)) and \
        len(list(x.glob('rains*')))>0]

# Me armo los inputs y outputs para el CDO

inputs = ' '.join([str(next(x.glob('rains*'))) for x in dir_physics])
output = str(dir_ensemble.joinpath('rains_every03_alltimes_wrfout_ensemble.nc'))

# CDO

CDO = cdo.Cdo()

CDO.ensmean(input=inputs, output=output)

