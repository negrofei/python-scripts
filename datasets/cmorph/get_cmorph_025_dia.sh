# Script para bajar y recortar CMORPH
echo Script para bajar y recortar CMORPH
read -p 'Elegí el año y mes (AAAAMM) ' AAAAMM
read -p 'Día que querés? ' DD
read -p 'Querés en 0.25 (1) o en 8km (2)? ' deg

MM=${AAAAMM:(-2)}
AAAA=${AAAAMM:0:4}

if [ $deg -eq 1 ]
then
    ruta='https://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V1.0/CRT/0.25deg-3HLY/'$AAAA'/'$AAAAMM'/CMORPH_V1.0_ADJ_0.25deg-3HLY_'$AAAAMM$DD'.bz2'
else
    ruta='https://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V1.0/CRT/8km-30min/'$AAAA'/CMORPH_V1.0_ADJ_8km-30min_'$AAAAMM'.tar'
    echo '### Mirá que para 8km está incompleto el script ###'
fi

MONTHS=(jan feb mar apr may jun jul aug sep oct nov dec)

wget $ruta
bzip2 -d *bz2

echo '### Hago el ctl ###'
rm cmorph.ctl || echo 'Creando el CTL file..'
echo 'DSET ^CMORPH_V1.0_ADJ_0.25deg-3HLY_'$AAAAMM$DD >> cmorph.ctl
echo 'OPTIONS little_endian
UNDEF  -999.0
TITLE  Precipitation estimates
XDEF 1440 LINEAR    0.125  0.25
YDEF  480 LINEAR  -59.875  0.25
ZDEF   01 LEVELS 1
TDEF 999999 LINEAR  00z'$DD${MONTHS[$MM-1]}$AAAA' 3hr
VARS 1
cmorph   1   99 RAW CMORPH integrated satellite precipitation estimates [mm/3hr]
ENDVARS
***  precipitation for the first time step every day (00Z) is the
***  accumulation for 00:00Z-02:59Z' >> cmorph.ctl


echo '### Convierto en NetCDF4 ###'
file='CMORPH_V1.0_ADJ_0.25deg-3HLY_'$AAAAMM$DD
cdo -f nc import_binary cmorph.ctl $file.nc

echo '### Recorto las lats lons ###'
cdo sellonlatbox,-70,-40.2,-40,-13.2 $file.nc ${file}_CSAM.nc
rm $file.nc
rm CMORPH_V1.0_ADJ_0.25deg-3HLY_$AAAAMM$DD

echo '### Listo papa ###'
