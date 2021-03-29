"""
 ### casos_mios ###

Calculo el FSS cada 3 horas a lo largo de cada evento. 
Esto se hace para cada miembro del ensamble y cada configuración 
Basicamente hago el calculo del FSS para todos los eventos, todas las configs y todas las simulacione


"""
# Paquetes
from netCDF4 import Dataset
import numpy as np
import pickle
import os
import sys
import time
import datetime
import Docpy
from Docpy.functions import printer
from tqdm import tqdm 
#####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Diccionarios y rutas
printer('Configuracion')

# inputs
acum = 3

# listas

datos_name = ['CMORPH','MSWEP','IMERG']

casos = [
         #'caso2005-03',
         'caso2015-11',
         #'caso2016-10',
         ]

configs = ['config_inicial',
           'igual1_con-sha',
           #'pbl_MYNN',
           #'shallow-on_MYNN_grell',
           #'shallow-on_MYNN_grell_MP',
           ]

sims = {
        #'2.4km_anidado12':'2',
        #'4km_anidado12':'2',
        #'4km_anidado12_B':'2',
        #'4km_anidado20_B':'2',
        #'4km_B':'1',
        '4km_B_15S47W':'1',
        }

# rutas
ruta_datos = os.path.abspath('./interpolated_files/')
try:
    os.mkdir('./pickles')
except:
    printer('Parece que ya existe, mono')
ruta_pickles = os.path.abspath('./backup_pickles')

ini = time.time()
for caso in tqdm(casos, desc='casos'):
    printer(caso)
    caso = caso[4:8]+caso[9:]
    casou = 'caso'+caso[:4]+'-'+caso[-2:]
    cmorph = Docpy.Dato(os.path.join(ruta_datos,'CMORPH_V1.0_ADJ_8km-30min_'+caso+'_CSAM_boxed.nc'))

    mswep = Docpy.Dato(os.path.join(ruta_datos,'mswep_'+caso+'_boxed.nc'))

    imerg = Docpy.Dato(os.path.join(ruta_datos,'3B-HHR.MS.MRG.3IMERG.'+caso+'_boxed.nc'))
    
    datos = [cmorph, mswep, imerg]

    common_times_datos = Docpy.functions.common_times([Docpy.functions.acum_time(dato, acum) for dato in datos])

    rutas = []
    for config in configs:
        for sim in sims:
            ruta = os.path.join(ruta_datos,casou,config,sim,'wrfout_arw_d0'+sims[sim]+'_'+casou[4:]+'_boxed.nc')
            if os.path.isfile(ruta):
                rutas.append(ruta)
                print(ruta)
    mascara = np.sum([Docpy.precip.calc_precip_acum(datos[i], acum)[0,:,:].mask for i in range(len(datos))], axis=0)>=1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Cargo las cosas
    printer('Cargo los datos y el modelo')
    for dato in datos:
        print(dato.name)
        data_precip_tot = Docpy.precip.calc_precip_acum(dato, acum=acum)
        acum_time_dato = Docpy.functions.acum_time(dato, acum)
        for ruta in rutas:
            print(ruta)
            model_class = Docpy.WRF(ruta)
            model_class.sim = ruta.split('/')[-2].replace('_','-')
            model_class.config = ruta.split('/')[-3]
            model = model_class.data_netcdf
            #%% Hago coincidir los tiempos en los que quiero calcular el FSS
            ### y selecciono la precip en esos tiempos
            acum_time_model = Docpy.functions.acum_time(model_class, acum)
            common_times = Docpy.functions.common_times([common_times_datos, acum_time_model])
            print(common_times)
            for tiempo in common_times:
                print(tiempo)
                data_precip = np.ma.masked_array(data_precip_tot[acum_time_dato.index(tiempo)], mascara)
                model_precip = np.ma.masked_array(Docpy.precip.calc_precip_acum(model_class, acum)[acum_time_model.index(tiempo)], mascara)
                
                ### ACA EMPIEZA EL FSS
                ny = model_precip.shape[0]
                nx = model_precip.shape[1]
                
                IM_WRF = np.zeros([ny,nx])
                IO_OBS = np.zeros([ny,nx])
                percentiles = [75,90,95,99]
                vector_n = np.array(np.arange(1,2*max(ny,nx)+1,2))
                M = np.zeros([int(np.max(vector_n))+1,ny,nx])
                O = np.zeros([int(np.max(vector_n))+1,ny,nx])
                MSE = np.zeros([int(np.max(vector_n))+1])
                MSE_ref = np.zeros([int(np.max(vector_n))+1])
                FSS = np.zeros([int(np.max(vector_n))+1])

#%% CALCULO DE PERCENTILES CERCANOS AL MÃ<81>XIMO

                for p in percentiles:
                    fss_name = '_'.join(['FSS',
                        'acum'+'{:02d}'.format(acum),
                        caso,
                        dato.name,
                        model_class.sim+'-'+model_class.config,
                        'P'+str(p),
                        tiempo])
                    if os.path.isfile(os.path.join(ruta_pickles,fss_name)):
                        print(fss_name,'ya existe!')
                    else:

                        per_model = np.nanpercentile(np.ma.array(model_precip).compressed(), p)
                        per_data = np.nanpercentile(np.ma.array(data_precip).compressed(), p)
                        printer(tiempo,char='-')
                        print('Percentil',p,
                              ' modelo',model_class.name, model_class.sim, model_class.config,
                              ':',round(per_model,2))
                        print('Percentil',p,' datos',dato.name,':',per_data)
                        
                        #%%%  CALCULO EL IO E IM PARA CADA TIEMPO
                        printer('#')
                        print(dato.name,'\t',model_class.name, '\t', model_class.sim,'\t',
                              model_class.dominio)
                        printer('#')
                        for j in range(ny):
                            for i in range(nx):
                                if model_precip[j,i]>=per_model:
                                    IM_WRF[j,i]=1
                                else:
                                    IM_WRF[j,i]=0

                                if data_precip[j,i]>=per_data:
                                    IO_OBS[j,i]=1
                                else:
                                    IO_OBS[j,i]=0
                        #%%% AHORA CALCULO LAS FRACCIONES SUMADAS
                        sumM=np.zeros([ny,nx])
                        sumO=np.zeros([ny,nx])
                        sumM[0,0]=IM_WRF[0,0] # seteo el primer valor como el IM en 0,0
                        sumO[0,0]=IO_OBS[0,0]
                        for j in range(1,ny): # genero la primera columna 
                            sumM[j,0]=sumM[j-1,0]+IM_WRF[j,0]
                            sumO[j,0]=sumO[j-1,0]+sumO[j,0]
                        for i in range(1,nx): #genero la primera fila
                            sumM[0,i]=sumM[0,i-1]+IM_WRF[0,i]
                            sumO[0,i]=sumO[0,i-1]+sumO[0,i]
                            for j in range(1,ny):
                                sumM[j,i]=sumM[j-1,i]+sumM[j,i-1]-sumM[j-1,i-1]+IM_WRF[j,i]
                                sumO[j,i]=sumO[j-1,i]+sumO[j,i-1]-sumO[j-1,i-1]+IO_OBS[j,i]
    #%% AHORA CALCULO EL FSS CON LAS FRACCIONES SUMADAS 
    ### LA FRACCIÃ<93>N REAL (M Y O) SE OBTIENEN AGARRANDO LAS FRACCIONES SUMADAS Y TOMANDO EL VALOR DE LA ESQUINA 
    ### DEL 'BARRIO', SE LE RESTA LAS ESQUINAS DE AFUERA Y SE LE SUMA EL VALOR DE LA ESQUINA OPUESTA DE AFUERA.
    ### VER EL TRABAJO DE JONES
                        for n in vector_n: #para cada n del barrio
                            for j in range(ny):
                                for i in range(nx): #para cada punto de reticula
                                    # PONGO LA CONDICION DE QUE NO SE ESCAPE DEL DOMINIO
                                    if i+(n-1)/2>=nx-1:
                                        imasn=int(nx-1)
                                    else:
                                        imasn=int(i+(n-1)/2)
                                    if i-(n-1)/2<=0:
                                        imenosn=0
                                    else:
                                        imenosn=int(i-(n-1)/2-1)
                                    if j+(n-1)/2>=ny-1:
                                        jmasn=int(ny-1)
                                    else:
                                        jmasn=int(j+(n-1)/2)
                                    if j-(n-1)/2<=0:
                                        jmenosn=0
                                    else:
                                        jmenosn=int(j-(n-1)/2-1)
                                    A=sumM[jmasn,imasn]
                                    a=sumO[jmasn,imasn]
                                    B=sumM[jmasn,imenosn]
                                    b=sumO[jmasn,imenosn]
                                    C=sumM[jmenosn,imasn]
                                    c=sumO[jmenosn,imasn]
                                    D=sumM[jmenosn,imenosn]
                                    if (i-(n-1)/2<=imenosn<=i+(n-1)/2) and (j-(n-1)/2<=jmasn<=j+(n-1)/2):
                                        B=D=b=d=0
                                    if (i-(n-1)/2<=imasn<=i+(n-1)/2) and (j-(n-1)/2<=jmenosn<=j+(n-1)/2):
                                        C=D=c=d=0
                                    if (i-(n-1)/2<=imenosn<=i+(n-1)/2) and (j-(n-1)/2<=jmenosn<=j+(n-1)/2):
                                        D=d=0
                                    # Los primeros 4 if's corresponden a ajustar la condicion de borde cuando la caja sobrepasa
                                    # el dominio. Los 3 if's del final son para que reste lo que tiene que restar.
                                    # BÃ¡sicamente, si algun punto que tomo en cuenta para la sumatoria estÃ¡ adentro de la caja,
                                    # lo vuelo y tmb vuelo el D, que estÃ¡ para ajustar nomÃ¡s.                
                                    M[int(n),j,i]=(A-B-C+D)/(n**2)
                                    O[int(n),j,i]=(a-b-c+d)/(n**2)
    #%%%  AHORA CALCULO TODO LO OTRO

                        for n in vector_n:
                            MSE[int(n)]=1/(ny*nx)*np.sum((M[int(n),:,:]-O[int(n),:,:])**2)
                            MSE_ref[int(n)]=1/(ny*nx)*(np.sum(M[int(n),:,:]**2+O[int(n),:,:]**2))

                        # FINALMENTE EL FSS
                        for n in vector_n:
                            if MSE_ref[int(n)]==0:
                                FSS[int(n)]=0
                            else:
                                FSS[int(n)]=1-MSE[int(n)]/MSE_ref[int(n)]

    #%% GUARDO LOS FSS y las IO e IM   
                        fss_name = '_'.join(['FSS',
                                             'acum'+'{:02d}'.format(acum),
                                             caso,
                                             dato.name,
                                             model_class.sim+'-'+model_class.config,
                                             'P'+str(p),
                                             tiempo])
                        with open(os.path.join(ruta_pickles,fss_name),'wb') as output:
                            pickle.dump(FSS,output,pickle.HIGHEST_PROTOCOL)

                        io_name = '_'.join(['IO',
                                            'acum'+'{:02d}'.format(acum),
                                            caso,
                                            dato.name,
                                            model_class.sim+'-'+model_class.config,
                                            'P'+str(p),
                                            tiempo])
                        with open(os.path.join(ruta_pickles,io_name),'wb') as output:
                            pickle.dump(IO_OBS,output,pickle.HIGHEST_PROTOCOL)

                        im_name = '_'.join(['IM',
                                            'acum'+'{:02d}'.format(acum),
                                            caso,
                                            dato.name,	
                                            model_class.sim+'-'+model_class.config,
                                            'P'+str(p),
                                            tiempo])
                        with open(os.path.join(ruta_pickles,im_name),'wb') as output:
                            pickle.dump(IM_WRF,output,pickle.HIGHEST_PROTOCOL)
        dato_time = time.time()
        printer('#')
        printer('-')
        printer(str(datetime.timedelta(seconds=int(dato_time-ini))))
    caso_time = time.time()
    printer('#')
    printer('-')
    printer(str(datetime.timedelta(seconds=int(caso_time-ini))))
fin = time.time()
print(str(datetime.timedelta(seconds=int(fin-ini))))

