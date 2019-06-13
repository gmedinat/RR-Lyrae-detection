#!/usr/bin/python2.7
 
import numpy as np#para cargar numpy
from pylab import *#cargar mtplotlib
import os#cargar herramientas de archivos
import re # use regular patterns
import sys, getopt # system commands
import string # string functions4
import math
#from scipy import linalg as scipylinalg
#from scipy import stats
#from scipy import ndimage
#from scipy import signal as scipysignal
from pylab import *
from numpy import nan
import random
import subprocess
#import lomb
import pyfits as fits
#from astroML.time_series import lomb_scargle
import glob
from astropy.stats import sigma_clip
import pandas as pd

# ESTE PROGRAMA VE LAS CURVAS DE LUZ EN ESTRELLAS GEN, Y FILTRA LAS QUE TENGAN COLORES G-R ENTRE -0.45 Y 0.45 
# las magnitudes en g las saca de los catalogos completos que estan en Data/fields_with_uncalibrated_mags


#indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s"%campo 
#files = os.listdir("/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s"%campo)

#print '\nField %s *********************************'%campo

# Field1_13_57994.14318666_r.csv
files_0 = glob.glob("/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J/*.csv")
#files_0 = glog.glob("Field_*/TimeSeries_r/Filtered_All_Periods_P4J/*.csv")

files_0.sort()

files = files_0
#files = files_0[37500:]

#print files
largo = len(files)

print largo
#sys.exit(0)

i=0
for f in files:

	print '\n'	
	
	#if f == '/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J/obj_306.2062306_-40.0276269_1_N1_r.csv':
	#	sys.exit(0)

	lc_df = pd.read_csv(f)
	print i+1, '/', largo, '\t', f

	#obj_311.0413557_-35.9696288_7_S6_r.csv

	# debemos extraer ra, dec, field, ccd: 4 elementos
	extracted_element = 0
	index = -7
	prev_index = -6
	while extracted_element <=4:
		if f[index] == '_':
			if extracted_element == 0:
				ccd = f[index+1:prev_index]
			if extracted_element == 1:
				field = f[index+1:prev_index]
			if extracted_element == 2:
				dec = f[index+1:prev_index]
			if extracted_element == 3:
				ra = f[index+1:prev_index]
				break
	
			prev_index = index
			extracted_element = extracted_element+1

		index = index-1
	
	ccdnum_useful = lc_df.iloc[0]['ccdnum']
	catalog_file = '/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Fields_with_uncalibrated_magnitudes/PerCCD/RR%s_ccdnum%s_g.csv'%(field,ccdnum_useful)
	cat_df = pd.read_csv(catalog_file)

	#cat_df = cat_df[['id', 'RA_deg', 'DEC_deg', 'coord_ra','coord_dec', 'object', 'Field', 'visit', 'ccdnum', 'MJDObs', 'Obs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1
	#'mag_rel_cal'


	# ver cuantas epocas hay (x)

	#unicos = cat_df.MJDObs.unique()
	#print unicos

	#n_epocas = len(unicos)

	# elegir los x mejores match

	# siesque es el mismo objeto ('object'), calcular promedio de magnitud

	# separar por epoca?

	# match

	ventana_match = 2./60./60. #  2 arcsec

	#df_sort = cat_df.ix[(((cat_df['RA_deg']-float(ra))**2.+(cat_df['DEC_deg']-float(dec))**2.)**0.5).abs().argsort()[:n_epocas]]
	#print df_sort

	cat_df['distancia_deg'] = (((cat_df['RA_deg']-float(ra))**2.+(cat_df['DEC_deg']-float(dec))**2.)**0.5)
	idxmin = cat_df['distancia_deg'].idxmin()
	
	#print df_obj_id

	#sys.exit(0)

	if cat_df.iloc[idxmin]['distancia_deg'] <= ventana_match:

		object_id = cat_df.iloc[idxmin]['object'] # == float(CCDnums[cont_CCD])
		df_obj_id = cat_df.loc[cat_df['object'] == object_id]

		print '\tlight curve source (ra, dec, field, ccd)\t    ||||||| \t catalog_g (ra, dec)\n\t', ra, dec, field, ccd, '\t \t \t \t', cat_df.iloc[idxmin]['RA_deg'], cat_df.iloc[idxmin]['DEC_deg']

		mag_g = df_obj_id['mag_rel_cal'].mean()
		mag_r = lc_df['mag_rel_cal'].mean()
		N_g = len(df_obj_id)
		N_r = len(lc_df)

		lc_df['object_g'] = str(object_id)
		lc_df['mean_mag_rel_cal_r'] = mag_r
		lc_df['mean_mag_rel_cal_g'] = mag_g
		lc_df['color'] = mag_g - mag_r
		lc_df['N_r'] = N_r
		lc_df['N_g'] = N_g
		lc_df['g_r_distance'] = cat_df.iloc[idxmin]['distancia_deg'] 

		if len(df_obj_id) >= 2:

			#mag_g = df_obj_id['mag_rel_cal'].mean()

			print 'mean_mag_g: ', mag_g, '\t object id: ', int(object_id), '\t Number of observations: ', len(df_obj_id), '\t Distance_deg: ', cat_df.iloc[idxmin]['distancia_deg']

			#mag_r = lc_df.iloc[0]['mag_rel_cal']
			#mag_r = lc_df['mag_rel_cal'].mean()
			print 'mean_mag_r: ', mag_r, '\t\t g-r: ', mag_g-mag_r

			if (mag_g-mag_r) >= -0.45 and (mag_g-mag_r) <= 1.0:

				print 'Color-cut selected!'
				#lc_df['mean_mag_rel_cal_r'] = mag_r
				#lc_df['mean_mag_rel_cal_r_err'] = 
				#lc_df['mean_mag_rel_cal_g'] = mag_g
				#lc_df['mean_mag_rel_cal_g_err'] = 
		
				lc_df.to_csv( f.replace('Estrellas_Gen_All_Periods_P4J','Estrellas_Gen_All_Periods_P4J/Color_cut_test') )
				#os.system( 'cp %s %s'%(f.replace('.csv','.png').replace('obj', 'starpngGLS') , f.replace('Estrellas_Gen_All_Periods_P4J','Estrellas_Gen_All_Periods_P4J/Color_cut').replace('.csv','.png').replace('obj', 'starpngGLS') ) )
				#break

				#if not os.path.exists(f.replace('.csv','.png')):
				#	os.system( 'cp %s %s'%(f.replace('Estrellas_Gen_All_Periods_P4J','Estrellas_Gen_All_Periods_P4J/Selected_first').replace('.csv','.png') , f.replace('Estrellas_Gen_All_Periods_P4J','Estrellas_Gen_All_Periods_P4J/Color_cut').replace('.csv','.png') ) )
	
				#else:
				#	os.system( 'cp %s %s'%(f.replace('.csv','.png') , f.replace('Estrellas_Gen_All_Periods_P4J','Estrellas_Gen_All_Periods_P4J/Color_cut').replace('.csv','.png') ) )
		elif len(df_obj_id) < 2:
			print 'ALERTA!! Best match has a unique object id. \n'
                        lc_df.to_csv( f )
		del df_obj_id		

	elif cat_df.iloc[idxmin]['distancia_deg'] > ventana_match:
		print 'ALERTA!! Match too far away. Best match distance: ', cat_df.iloc[idxmin]['distancia_deg'], '\n'


	#unicos_obj = df_sort.object.unique()
	#n_obj = len(unicos_obj)
	#if n_obj > 1:
	#	print '\n\tALERTA!!\t numero de objetos que minimizan distancia: ', n_obj

	del lc_df
	del cat_df

	#del df_sort
	
	#if i%500 == 0:
                # inform which was the last i%1000 == 0 that was examined. np arange 0,2,1 is only to put some content
      #          np.savetxt('checkpoints/checkpoint_filtro_color_%i_%s.check'%(i,largo), np.arange(0.0,2.0,1.0))
      #          if i != 0:
      #                  i_anterior = i-500
      #                  os.system( 'rm checkpoints/checkpoint_filtro_color_%i_%s.check'%(i_anterior, largo) )


	i = i+1


	
