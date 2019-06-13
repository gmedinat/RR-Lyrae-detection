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

#table = Table.read('testRR1_all_g.fits', format='fits')
#df = table.to_pandas()

#print df.columns()
#print df['id'] # la columna id
#print df.head(10) #las primeras 10 filas

#print df['id'] + df['MJDObs']

#print astropy.__version__


#df['nombrenuevo']=df['id'] + df['MJDObs']

#print df['MJDObs'].sort_values() # sort_values
#df.sort_values(by='MJDObs')

#diff_epochs = df.MJDObs.unique()
#print diff_epochs, len(diff_epochs)

#epoca = 1
#i=0
#while i < len(diff_epochs):
#	sub_tabla = df.loc[df['MJDObs'] == diff_epochs[i]]
#	df.to_csv('hola_%s.csv'%epoca)

#	print epoca

#	epoca = epoca+1
#	i=i+1

start = timeit.default_timer()

indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data" 

files = os.listdir("/run/media/gmedina/TOSHIBA EXT/DECam2017/Data")

#print files
#newfiles = files[47:]
#print "\n", newfiles, '\n'
#newfiles = files

#sys.exit(0)

modo_banda_gyr = False

if modo_banda_gyr == True: # para hacer ambas bandas al mismo tiempo

	file_n = 1
	band = 'g'
	count = 1

	file_n = 15
	band = 'r'
	count = 30

if modo_banda_gyr == False:
	file_n = 6
	#band = 'g'
	band = 'r'
	count = 6

print '\n'

#for f in newfiles:
while count <= 16 :

	#file names: RR9_all_g.fits, RR14_all_g.fits
	#if band_n == 0:
	#	band = 'g'
	#else:
	#	band = 'r'
	archivo = 'RR%i_all_%s.csv'%(file_n,band)

	# testRR1_all_g.dat
	#f = 'testRR1_all_r.fits'
	# RR10_all_withmd_g.fits

	#if f[len(f)-4:]=='fits':

	#	if f[3] == '_':
	#		field_n = f[2]
	#	else:
	#		field_n = f[2:4]
	#	band    =   f[len(f)-6]

	#	if file_n%2 != 0:
	#		banda = 'g'
	#	else:
	#		banda = 'r'

	#	archivo = f.replace(str(field_n), str(file_n))
	#	archivo = archivo.replace(band+'.fits', banda+'.fits')
		#outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s"%field_n
	outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s"%file_n

		#print f, '\t', field_n, band
	print archivo, '\t', file_n, band

	# Crear carpeta para este campo
	#if os.path.exists(outdir):
	#    os.system( 'rm %s/*'%(outdir) ) # eliminar los datos anteriores del out
	#if not os.path.exists(outdir):
	#    os.makedirs(outdir)

	#table = Table.read(indir+'/'+archivo, format='fits')
	#df = table.to_pandas()

	#df = pd.read_csv(indir+'/RR%s/'%file_n+archivo)
	cols = ['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']
	df = pd.read_csv(indir+'/Field_%s/'%file_n+archivo, usecols=cols)

	coord_ra  = df['coord_ra']
	coord_dec = df['coord_dec']
	df['RA_deg']  = coord_ra*180./math.pi
	df['DEC_deg'] = coord_dec*180./math.pi

	df = df[['id','coord_ra','coord_dec', 'RA_deg', 'DEC_deg', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1

	#df.sort_values(by='MJDObs')
	diff_epochs = df.MJDObs.unique()
	diff_epochs = np.sort(diff_epochs)
	print "\t epocas: ", diff_epochs, "\t epocas totales: ", len(diff_epochs)

	#break # las epocas estan ordenadas?

	epoca = 1
	i=0
	while i < len(diff_epochs):
		sub_tabla = df.loc[df['MJDObs'] == diff_epochs[i]]
		#sub_tabla = sub_tabla[['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1
		#sub_tabla.to_csv(outdir+"/Field%s_%s_%s_%s.csv"%(field_n,epoca, diff_epochs[i], band))
		sub_tabla.to_csv(outdir+"/Field%s_%s_%s_%s.csv"%(file_n,epoca, diff_epochs[i], band))

		print '\t \t \t', epoca

		epoca = epoca+1
		i=i+1

	del df
	del sub_tabla

	if modo_banda_gyr == True:
		if count%2==0:	
			file_n = file_n+1
		if band == 'g':
			band = 'r'
		elif band == 'r':
			band = 'g'
	if modo_banda_gyr == False:
		file_n = file_n+1
	count = count+1


stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))


#df['object']
#print list(df.columns.values)
#print list(df)


#fixed_df = pd.read_csv('testRR1_all_g.csv')
#print fixed_df[:3]

#print fixed_df['object']
#fixed_df['object'].plot()

#from astropy.table import Table

#table = Table.read('testRR1_all_g.fits')
#pandas_df = table.to_pandas()

#fixed_df = pd.read_csv('testRR1_all_g.csv', sep=';', encoding='latin1', parse_dates=['Date'], dayfirst=True, index_col='Date')
#print fixed_df[:3]

#img_path = "/run/media/gmedina/TOSHIBA\ EXT/DECam2017/Data"

#fits_image_filename = fits.util.get_testdata_filepath('%s/testRR1_all_g.fits'%img_path)
#fits_image_filename = fits.util.get_testdata_filepath('./testRR1_all_g.fits')
#hdul = fits.open('./testRR1_all_g.fits')

#print hdul.info()
#print hdul.data.shape

#print hdul[0].data.shape
