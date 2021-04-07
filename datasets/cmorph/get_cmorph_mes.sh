# Script para bajar y recortar CMORPH en la caja HR3
#
# Devuelve archivos CMORPH_V1.0_ADJ_<horizontal_resolution>_<AAAAMM>_CSAM1.nc

echo Script para bajar y recortar CMORPH
read -p 'Elegí el año y mes (AAAAMM) ' AAAAMM
read -p 'Querés en 0.25 (1) o en 8km (2)? ' deg

MM=${AAAAMM:(-2)}
AAAA=${AAAAMM:0:4}


MONTHS=(jan feb mar apr may jun jul aug sep oct nov dec)
dias=(31 28 31 30 31 30 31 31 30 31 30 31)

dia=${dias[(10#$MM-1)]}

if [ $deg -eq 1 ]
then

    for DD in $(seq -w 01 $dia)
    do  
        echo $DD    
        echo '### Hago el ctl ###'
        ruta='https://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V1.0/CRT/0.25deg-3HLY/'$AAAA'/'$AAAAMM'/CMORPH_V1.0_ADJ_0.25deg-3HLY_'$AAAAMM$DD'.bz2'
        wget $ruta
        bzip2 -d -v *bz2
        rm cmorph.ctl || echo 'Creando el CTL file..'
        echo 'DSET ^CMORPH_V1.0_ADJ_0.25deg-3HLY_'$AAAAMM$DD >> cmorph.ctl
        echo 'OPTIONS little_endian
        UNDEF  -999.0
        TITLE  Precipitation estimates
        XDEF 1440 LINEAR    0.125  0.25
        YDEF  480 LINEAR  -59.875  0.25
        ZDEF   01 LEVELS 1
        TDEF 999999 LINEAR  00z'$DD${MONTHS[10#$MM-1]}$AAAA' 3hr
        VARS 1
        cmorph   1   99 RAW CMORPH integrated satellite precipitation estimates [mm/3hr]
        ENDVARS
        ***  precipitation for the first time step every day (00Z) is the
        ***  accumulation for 00:00Z-02:59Z' >> cmorph.ctl
           
        echo '### Convierto en NetCDF4 ###'
        file='CMORPH_V1.0_ADJ_0.25deg-3HLY_'$AAAAMM$DD
        cdo -f nc import_binary cmorph.ctl $file.nc
        cdo seltimestep,1,2,3,4,5,6,7,8 $file.nc ${file}2.nc
        echo '### Recorto las lats lons ###'
        cdo sellonlatbox,-76,-48,-36,-16 ${file}2.nc ${file}_CSAM.nc
        rm $file.nc
        rm ${file}2.nc
        rm CMORPH_V1.0_ADJ_0.25deg-3HLY_$AAAAMM$DD
    done
    echo '#### Concateno todo en un solo archivo ####'
    cdo cat CMORPH_V1.0_ADJ_0.25deg-3HLY_${AAAAMM}*_CSAM.nc CMORPH_V1.0_ADJ_0.25deg-3HLY_${AAAAMM}_CSAM1.nc
    rm *CSAM.nc
    echo '### Listo papa ###'

elif [ $deg -eq 2 ]
then
    ruta='https://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V1.0/CRT/8km-30min/'$AAAA'/CMORPH_V1.0_ADJ_8km-30min_'$AAAAMM'.tar'
    if test -f "CMORPH_V1.0_ADJ_8km-30min_"${AAAAMM}".tar"
    then
        echo "CMORPH_V1.0_ADJ_8km-30min_"${AAAAMM}".tar existe"
    else
        echo "##### Bajo el archivo #####"
        wget $ruta
    fi
    tar -xvf CMORPH_V1.0_ADJ_8km-30min_${AAAAMM}.tar
    cd ${AAAAMM}
    for DD in $(seq -w 01 $dia)
    do
        bzip2 -d -v 'CMORPH_V1.0_ADJ_8km-30min_'${AAAAMM}${DD}*.bz2
    done
    rm *bz2

    echo '##### Hago los CTLs y paso los archivos a NETCDF #####'
    for DD in $(seq -w 01 $dia)
    do
        for HH in $(seq -w 00 23)
        do
            rm cmorph.ctl || echo 'Creando el CTL file...'
            echo 'DSET ^CMORPH_V1.0_ADJ_8km-30min_'$AAAAMM$DD$HH >> cmorph.ctl
            echo 'OPTIONS little_endian
UNDEF  -999.0
TITLE  Precipitation estimates
XDEF 4948 LINEAR   0.036378335 0.072756669
YDEF 1649 LINEAR -59.963614    0.072771377
ZDEF   01 LEVELS 1
TDEF 999999 LINEAR  '$HH'z'$DD${MONTHS[10#$MM-1]}$AAAA' 30mn
VARS 1
cmorph   1  99  hourly cmorph [ mm/hr ]
ENDVARS' >> cmorph.ctl
        
            echo '##### Convierto en NetCDF #####'
            file='CMORPH_V1.0_ADJ_8km-30min_'$AAAAMM$DD$HH
            cdo -f nc import_binary cmorph.ctl $file.nc
            cdo seltimestep,1,2 $file.nc ${file}2.nc
            echo '##### Recorto las lats lons #####'
            cdo sellonlatbox,-76,-48,-36,-16 ${file}2.nc ${file}_CSAM.nc
            rm ${file}.nc
            rm ${file}2.nc
            rm CMORPH_V1.0_ADJ_8km-30min_$AAAAMM$DD$HH
        done
    done
    echo '##### Concateno todo en un solo archivo #####'
    cdo cat CMORPH_V1.0_ADJ_8km-30min_${AAAAMM}*_CSAM.nc CMORPH_V1.0_ADJ_8km-30min_${AAAAMM}_CSAM1.nc
    rm *CSAM.nc
    mv *CSAM1.nc ..
    cd ..
    rm -r ${AAAAMM}
    rm "*"${AAAAMM}".tar" 
    echo '### Listo papa ###'
fi
