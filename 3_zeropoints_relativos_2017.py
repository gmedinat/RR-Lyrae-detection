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
import timeit
#para el campo X, usa stilts para hacer match entre la epoca de ref y las demas. Al match le calcula las magnitudes y las diferencias en magnitudes. Para cada epoca, guarda un promedio, que luego se guarda en un archivo para cada campo, en el formato:
# epoch mag_diff(mag_ref-mag_epo)
# 01  0.1
# 03  -0.2
# 04  0.009


#__all__ = ['sigma_clip', 'binom_conf_interval', 'binned_binom_proportion']
#def sigma_clip(data, sig=3, iters=1, cenfunc=np.median, varfunc=np.var, maout=False):
#	data = np.array(data, copy=False)
#	oldshape = data.shape
#	data = data.ravel()
#	mask = np.ones(data.size, bool)
#	if iters is None:
#		i = -1
#		lastrej = sum(mask) + 1
#		while(sum(mask) != lastrej):
#			i += 1
#			lastrej = sum(mask)
#			do = data - cenfunc(data[mask])
#			mask = do * do <= varfunc(data[mask]) * sig ** 2
#		iters = i + 1
#	#TODO: ?print iters to the log if iters was None?
#	else:
#		for i in range(iters):
#			do = data - cenfunc(data[mask])
#			mask = do * do <= varfunc(data[mask]) * sig ** 2
#	if maout:
#		return np.ma.MaskedArray(data, ~mask, copy=maout != 'inplace')
#	else:
#		return data[mask], mask.reshape(oldshape)

start = timeit.default_timer()

#filtro = 'g'
filtro = 'r'

if filtro == 'r':
	epoca_ref = '8'
if filtro == 'g':
	epoca_ref = '1'

Matriz = None

campo = 1
while campo <= 16:

	indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s"%campo 
	files = os.listdir("/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s"%campo)

	print '\nField %s *********************************'%campo

	# Field1_13_57994.14318666_r.csv
	catalogo_ref = glob.glob("/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s/Field%s_%s_*_%s.csv"%(campo,campo, epoca_ref,filtro))[0]

	#catalogo_ref = catalogo_ref.replace(indir+'/', '')
	#df_ref = pd.read_csv(catalogo_ref)

	outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s/zp_files"%campo
	# Crear carpeta para este campo
	#if os.path.exists(outdir):
	#    os.system( 'rm %s/*'%(outdir) ) # eliminar los datos anteriores del out
	if not os.path.exists(outdir):
	    os.makedirs(outdir)

	Matriz = None

	first_time = 100 # solo un indicador para despues
	i = 1
	for f in files:

		if f[len(f)-5:]=='%s.csv'%filtro and f[len(f)-9:len(f)-6] != 'all':
			
			if f[6]=='_' and f[8]=='_':
				epoca = f[7]
				#busqueda_guion
				j = 9
				while f[j]!='_':
					j=j+1
				MJD = f[9:j]
			if f[7]=='_' and f[9]=='_':
				epoca = f[8]
				#busqueda_guion
				j = 10
				while f[j]!='_':
					j=j+1
				MJD = f[10:j]
			if f[6]=='_' and f[9]=='_':
				epoca = f[7:9]
				#busqueda_guion
				j = 10
				while f[j]!='_':
					j=j+1
				MJD = f[10:j]
			if f[7]=='_' and f[10]=='_':
				epoca = f[8:10]
				#busqueda_guion
				j = 11
				while f[j]!='_':
					j=j+1
				MJD = f[11:j]

			#table = Table.read(indir+'/'+f, format='fits')
			#df = table.to_pandas()
			#df = pd.read_csv(indir+'/'+f)

			#xx=df['base_SdssShape_psf_xx']
			#yy=df['base_SdssShape_psf_yy']
			#df['rad'] = 0.263*np.sqrt(xx**2 + yy**2)   # where 0.263''/pix is the plate scale of DECam
			#mean = df['rad'].mean()

			print '\t \t', f, '\t -------------- \t', catalogo_ref.replace(indir+'/', '')

			#os.system('java -jar stilts.jar tmatch2 matcher=2d in1=%s/%s ifmt1=ascii in2=%s/%s ifmt2=ascii out=%s/tabla.dat ofmt=ascii values1="X Y" values2="X Y" params=2'%(indir, referencia, indir, obser, indir) )

			if catalogo_ref.replace(indir+'/', '') == f:
				f=f
				#dat = np.array([[MJD, epoca, mean]])
				#if Matriz == None:
				#if i == 1:
				#	Matriz = dat
				#else:
				#	Matriz = np.vstack([Matriz, dat])
				MJD_ref = MJD
				if os.path.exists('/run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data/Field_%s/tabla.csv'%campo):
				    os.system('rm /run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data/Field_%s/tabla.csv'%campo)
				print '\t \t \t Done!'
			else:
				os.system('java -jar stilts.jar tskymatch2 in1=/run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data/Field_%s/%s ifmt1=csv in2=/run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data/Field_%s/%s ifmt2=csv out=/run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data/Field_%s/tabla.csv ofmt=csv ra1="RA_deg" dec1="DEC_deg" ra2="RA_deg" dec2="DEC_deg" error=1 join=1and2'%(campo, catalogo_ref.replace(indir+'/', ''), campo, f, campo) )
				#abrir tabla
				catalog_df = pd.read_csv(indir+'/tabla.csv')
				#catalog_df = pd.read_fwf('/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s/tabla.dat'%campo)         
				#print catalog_df.columns()

				catalog_df['a_1'] = catalog_df['ccdnum_1']
				catalog_df['b_1'] = catalog_df['ccdnum_1']
				catalog_df['k_1'] = catalog_df['ccdnum_1']
				
				catalog_df['a_2'] = catalog_df['ccdnum_2']
				catalog_df['b_2'] = catalog_df['ccdnum_2']
				catalog_df['k_2'] = catalog_df['ccdnum_2']

				# crear columnas a_1, k_1, a_2, k_2

				data_decam = np.loadtxt( '/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/DECAM_zps_%s.dat'%filtro , dtype='str')
				aas     = data_decam[:,2]
				aas_err = data_decam[:,3]
				bbs     = data_decam[:,4]
				bbs_err = data_decam[:,5]
				kks     = data_decam[:,6]
				kks_err = data_decam[:,7]

				cont = 1
				while cont <= 62:
					a = aas[cont-1]
					k = kks[cont-1]
					b = bbs[cont-1]

					a_err = aas_err[cont-1]
					b_err = bbs_err[cont-1]
					k_err = kks_err[cont-1]

					#donde ccdnum_1 sea "cont", que a_1 valga a, y k_1 valga k. 

					#catalog_df.a_1.fillna(df.Farheit, inplace=True)
					#catalog_df.a_1.fillna(a, inplace=True)
					#catalog_df.a_1.replace(to_replace=cont, a, inplace=True)
					#catalog_df.k_1.replace(to_replace=cont, k, inplace=True)
					catalog_df['a_1'].replace(cont, float(a), inplace=True)
					catalog_df['b_1'].replace(cont, float(b), inplace=True)
					catalog_df['k_1'].replace(cont, float(k), inplace=True)

					catalog_df['a_2'].replace(cont, float(a), inplace=True)
					catalog_df['b_2'].replace(cont, float(b), inplace=True)
					catalog_df['k_2'].replace(cont, float(k), inplace=True)

					cont = cont+1

				#print list(catalog_df.columns.values)
				#catalog_df.to_csv(indir+'/tabla.csv')				
				#sys.exit(0)

				#g = g_instr - a_g - b_g*( (g-r) - (g-r)0 ) - k_g * X
				
				#readfile
				#a = 
				#b = 
				#k = 
				
				X1 = catalog_df['AIRMASS_1']
				T1 = catalog_df['EXPTIME_1']
				a1 = catalog_df['a_1']
				k1 = catalog_df['k_1']

				X2 = catalog_df['AIRMASS_1']
				T2 = catalog_df['EXPTIME_2']
				a2 = catalog_df['a_2']
				k2 = catalog_df['k_2']
				

				f1 = catalog_df['base_PsfFlux_flux_1']
				catalog_df['mag_1'] = -2.5*np.log10(f1/T1) - a1 - k1*X1

				f2 = catalog_df['base_PsfFlux_flux_2']
				catalog_df['mag_2'] = -2.5*np.log10(f2/T2) - a2 - k2*X2

				catalog_df['mag_diff'] = catalog_df['mag_1'] - catalog_df['mag_2']

				catalog_df=catalog_df[['coord_ra_1','coord_dec_1', 'RA_deg_1', 'DEC_deg_1', 'base_PsfFlux_flux_1', 'base_PsfFlux_fluxSigma_1', 'base_PsfFlux_flux_2', 'base_PsfFlux_fluxSigma_2', 'mag_1', 'mag_2', 'mag_diff']]#1
				#print catalog_df.columns

				#print catalog_df['mag_diff']
				


				#sigma clipping
                        #datata, mask = sigma_clip(mag_diff, 2., None, mean) # returns: filtered data (array), mask (boolean array). 2.: std*2.    
				#magOBSclipp = magOBS[mask]
				#magREFclipp = magREF[mask]
				#diff_clipped = mag_diff[mask]
                        
				#diff = np.median(diff_clipped)
                        #calcular error. Basado en las magnitudes menores a 20.5?
				#mascara_errores = magOBSclipp < 20.5
				#sigma = np.std( diff_clipped[mascara_errores] )
				#sigmaN = sigma/math.sqrt(len(diff_clipped[mascara_errores]))

				#sacarles los nans
				#catalog = catalog[~np.isnan(catalog).any(axis=1)]

				filtered_data = sigma_clip(catalog_df['mag_diff'], sigma=3, iters=None)
				mask = ~filtered_data.mask
				#print filtered_data
				#print mask
				#print filtered_data[mask]
				#print df[mask]

				df_clip = catalog_df[mask]

				#dat = np.array([[MJD, epoca, mean]])
				#if Matriz == None:
				#if i == 1:
				#	Matriz = dat
				#else:
				#	Matriz = np.vstack([Matriz, dat])

				magOBS = catalog_df['mag_2']	
				mag_diff = catalog_df['mag_diff']

				magOBS_clip = df_clip['mag_2']	
				mag_diff_clip = df_clip['mag_diff']

				median = catalog_df['mag_diff'].median()
				median_clip = df_clip['mag_diff'].median()
				print median, median_clip

				sigma = df_clip['mag_diff'].std()
				sigmaN = sigma/math.sqrt(len(df_clip['mag_diff']))

				#guardar datos para hacer el grafico
				tabla_guardar = catalog_df[['mag_2','mag_diff']]
				tabla_guardar.to_csv(outdir+"/file_mag_diffs_%s_%s.csv"%(epoca, epoca_ref))

				tabla_guardar_clip = df_clip[['mag_2','mag_diff']]
				tabla_guardar_clip.to_csv(outdir+"/file_mag_diffs_%s_%s_clipped.csv"%(epoca, epoca_ref))

				#armar matriz de  pto zeros relativos para este campo
				dat = np.array([[MJD, epoca, median_clip, sigmaN]])
                        #if Matriz == None:
				if first_time == 100:
					Matriz = dat 					
				else:
					Matriz=np.vstack([Matriz, dat])

				#graficar magnitudes en Obs vs diferencia
				fig = plt.figure()
				ax = fig.add_subplot(111)
				ax.plot(magOBS, mag_diff, 'ro', label='all')
				ax.plot(magOBS_clip, mag_diff_clip, 'bo', label='clipped')
				#ax.plot(flujosREF*data[0], flujosOBS, 'bo', label='REF*aflux vs OBS')
				ax.plot(magOBS ,magOBS*0+median, 'r-')
				ax.plot(magOBS ,magOBS*0+median_clip, 'b-')
				ax.set_xlabel('$g_{OBS}$')
				ax.set_ylabel('$g_{REF}\ -\ g_{OBS}$')
				#ax.set_xlim(15, 24)
				#ax.set_ylim(-0.8, 0.8)
				ax.set_title('median: %f ||| median_clip: %f'%(median, median_clip) )
				plt.savefig( outdir+'/mag_diffs_%s_%s_%s.png'%(epoca, epoca_ref, filtro) )
				plt.close()
			
				os.system('rm /run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data/Field_%s/tabla.csv'%campo)

				print '\t \t \t Done!'

				del catalog_df
				del df_clip
				del tabla_guardar
				del tabla_guardar_clip
				first_time = 200	
				
			i = i+1
	
	dat_ref = np.array([[MJD_ref, epoca_ref, 0, 0]])
	Matriz=np.vstack([Matriz, dat_ref])	
	#ordenar Matriz segun la 1ra columna
	Matriz=Matriz[Matriz[:,0].argsort()]

	#agregar columna de MJDs
	#MJD = np.asarray(MJDs)
	
	#MJD = np.reshape( MJD, (len(MJD), 1) )


	#Matriz = np.hstack((MJD, Matriz))

	#guardar la matriz en la carpeta en un archivo de texto, de modo que matchPY2magPrevGen.py la pueda 	leer
	np.savetxt( outdir+"/zp_relativos_%s_%s.dat"%(campo, filtro), Matriz, fmt='%s',delimiter="   ", header='MJD \t Obs \t magRef-magObs \t error')

	campo = campo+1
	


stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))












