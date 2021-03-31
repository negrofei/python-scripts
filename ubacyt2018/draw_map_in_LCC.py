"""

Script for plotting a map of a simulation made in Lambert Conformal (LCC) projection.
This is based in the output files the Support team is developing to process WRF simulations

for running

$ python3 draw_map_in_LCC.py <path_to_file> 
"""

import sys                          # to read the file
import numpy as np                  # alltime numpy 
from netCDF4 import Dataset         # to read the netcdf file
import cartopy.crs as ccrs          # to get the projections
import cartopy.feature as cfeature  # to make the borderlines of countries
import matplotlib.pyplot as plt     # to plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Reading the file

filepath = sys.argv[1]
data = Dataset(filepath)

# Variables and dimensions
variables = [*data.variables]
variables.remove('Lambert_Conformal')       # remove projection info 
dimensions = [*data.dimensions]

selected_variables = [var for var in variables if var not in dimensions]

if len(selected_variables)==1:
    var = selected_variables[0]
else:
    var = input('Choose from the following:\n'+'\n'.join(selected_variables)+'\n')

# Select fields

dat = data.variables[var]
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
lev = data.variables['lev'][:]
x = data.variables['x'][:]*1000     # *1000 to convert to m instead of km
y = data.variables['y'][:]*1000

x2d, y2d = np.meshgrid(x, y)        # make meshgrid to allow transformation later

# Projection management

wrfproj = data.variables['Lambert_Conformal']       # load the projection of wrf simulation
lambert = ccrs.LambertConformal(                    # load the LCC for cartopy
        central_longitude=float(wrfproj.longitude_of_central_meridian),
        central_latitude=float(wrfproj.latitude_of_projection_origin),
        standard_parallels=(float(wrfproj.northern_parallel), float(wrfproj.southern_parallel)),
        cutoff=np.max(lat)+1
        )

# Transform from latlon to LCC

transform = ccrs.PlateCarree().transform_points(lambert, x2d, y2d)
x_lcc = transform[...,0]
y_lcc = transform[...,1]

# Plot
bordes = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        edgecolor='k',
        facecolor='none'
        )

while True:
    print('Set time and lev as in Grads:')
    times = input('Time: ')
    levs = input('Lev: ')
    
    t_idx = int(times.split(' ')[-1])
    lev_idx = list(lev).index( int(levs.split(' ')[-1]) )
    fig = plt.figure()
    axes = plt.axes(projection=lambert)
    mapa = axes.contourf(x_lcc, y_lcc, dat[t_idx, lev_idx,:],
            transform=ccrs.PlateCarree())

    # Config
    axes.add_feature(bordes)
    axes.coastlines('50m', color='k')
    
    plt.colorbar(mapa)
    plt.show()
    
    while True:
        cont = input('Continue? ([Y]/n) ')
        if cont.upper() == 'Y' or cont == '':
            no = False
            break
        elif cont == 'n':
            no = True
            break
        else:
            print('Not an answer')

    if no:
        break




