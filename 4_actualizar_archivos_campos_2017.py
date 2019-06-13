#!/usr/bin/python2.7
 
import numpy as np#para cargar numpy
from pylab import *#cargar mtplotlib
import os#cargar herramientas de archivos
import re # use regular patterns
import sys, getopt # system commands
import string # string functions4
import math
from pylab import *
from numpy import nan
import random
import subprocess
import csv 
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from astropy.io import fits 
import pandas as pd

from astropy.table import Table
import astropy
import timeit

start = timeit.default_timer()


#band = 'g'
band = 'r'

indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data" 

field = 15
while field <= 16 :


	#archivo = 'RR%i_all_%s.fits'%(field,band)
	archivo = 'RR%i_all_%s.csv'%(field,band)

	outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Fields_with_uncalibrated_magnitudes/"
	
	print archivo, '\t', field, band

	if not os.path.exists(outdir):
	    os.makedirs(outdir)

	#table = Table.read(indir+'/'+archivo, format='fits')
	#df = table.to_pandas()
	cols = ['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']
	df = pd.read_csv(indir+'/Field_%i/'%field+archivo, usecols=cols)
	#df = df[['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1


	# ------------ Calcular magnitudes segun formula que usea para los zp relativos (DECam zeropoints) ----------

	df['a'] = df['ccdnum']
	df['b'] = df['ccdnum']
	df['k'] = df['ccdnum']

	df['a_err'] = df['ccdnum']
	df['b_err'] = df['ccdnum']
	df['k_err'] = df['ccdnum']
				
	# crear columnas a, k
	data_decam = np.loadtxt( '/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/DECAM_zps_%s.dat'%band , dtype='str')
	aas     = data_decam[:,2]
	aas_err = data_decam[:,3]
	bbs     = data_decam[:,4]
	bbs_err = data_decam[:,5]
	kks     = data_decam[:,6]
	kks_err = data_decam[:,7]

	cont = 1
	while cont <= 62:
		a_val = aas[cont-1]
		k_val = kks[cont-1]
		b_val = bbs[cont-1]

		a_err = aas_err[cont-1]
		b_err = bbs_err[cont-1]
		k_err = kks_err[cont-1]

		#donde ccdnum sea "cont", que 'a' valga a, y 'k' valga k. 

		df['a'].replace(cont, float(a_val), inplace=True)
		df['b'].replace(cont, float(b_val), inplace=True)
		df['k'].replace(cont, float(k_val), inplace=True)

		df['a_err'].replace(cont, float(a_err), inplace=True)
		df['b_err'].replace(cont, float(b_err), inplace=True)
		df['k_err'].replace(cont, float(k_err), inplace=True)

		cont = cont+1

	X = df['AIRMASS']
	T = df['EXPTIME']
	a = df['a']
	k = df['k']
	aerr = df['a_err']
	berr = df['b_err']
	kerr = df['k_err']


	f = df['base_PsfFlux_flux']
	ef = df['base_PsfFlux_fluxSigma']	

	mag_name = 'mag_%s'%band
	df[mag_name] = -2.5*np.log10(f/T) - a - k*X # le pongo mag_"filtro"?


	print '\t\t Magnitudes instrumentales agregadas exitosamente!'

	# ------------- add zeropoint relativo --------------------
	print '\t\t\t abriendo zp_relativos_%s_%s.dat'%(field,band)
	data_zp_rel = np.loadtxt( '/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s/zp_files/zp_relativos_%s_%s.dat'%(field,field,band) , dtype='str')
	MJD_rel     	= data_zp_rel[:,0]
	Obs_rel 		= data_zp_rel[:,1]
	ref_min_epoch     = data_zp_rel[:,2]
	error 		= data_zp_rel[:,3]


	df['zp_rel'] = df['MJDObs']
	df['zp_rel_err'] = df['MJDObs']
	df['Obs'] = df['MJDObs']

	n_epochs = len(MJD_rel)
	count = 1
	while count <= n_epochs:
		MJD_val 	= MJD_rel[count-1]
		Obs_val 	= Obs_rel[count-1]
		ref_min_val = ref_min_epoch[count-1]
		err_val 	= error[count-1]

		#donde ccdnum sea "cont", que 'a' valga a, y 'k' valga k. 
		
		df['zp_rel'].replace(float(MJD_val), float(ref_min_val), inplace=True)
		df['zp_rel_err'].replace(float(MJD_val), float(err_val), inplace=True)
		df['Obs'].replace(float(MJD_val), float(Obs_val), inplace=True)

		count = count+1

	print '\t\t Zeropoint relativos agregados exitosamente!'

	mag    = df[mag_name]
	zp_rel = df['zp_rel']
	edeltaMag_zp = df['zp_rel_err']

	df['mag_rel_cal'] = mag + zp_rel
	df['mag_rel_cal_err'] = np.sqrt( (2.5*ef/f/np.log(10))**2 + aerr**2 + (X*kerr)**2 + (edeltaMag_zp)**2 )

	ra = df['coord_ra']
	dec = df['coord_dec']

	df['RA_deg']  = ra*180./math.pi
	df['DEC_deg'] = dec*180./math.pi

	df['Field']   = field

	df = df[['RA_deg', 'DEC_deg','coord_ra','coord_dec', 'object', 'Field', 'ccdnum', 'MJDObs', 'Obs', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma', 'a', 'a_err', 'b', 'b_err', 'k', 'k_err', mag_name, 'zp_rel', 'zp_rel_err', 'mag_rel_cal','mag_rel_cal_err']]#1

	df.to_csv(outdir+"/RR%s_%s.csv"%(field,band))
	
	del df

	field = field+1

 
stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))














