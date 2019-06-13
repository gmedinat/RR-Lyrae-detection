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
import corner
import scipy.optimize as op
import emcee
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

from astropy.io import fits 
import pandas as pd

from astropy.table import Table
import astropy 


# los catalogos en g ya estan juntados, pero estan en formato fits. Hay que cambiarlos a cvs.
# los catalogos en r hay que juntarlos.

#band = 'g'
band = 'r'


indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data" 

print '\n\n'

field = 1

while field <= 24 :
	outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/RR%i"%field 
	
	print "FIELD ", field, '------------------------------------------------------------------'

	if band == 'x': # solo convertir el formato
		# RR20_all_g.fits
		archivo = 'RR%i/RR%i_all_%s.fits'%(field,field,band)
		table = Table.read(indir+'/'+archivo, format='fits')
		stacked_df = table.to_pandas()

	if band == 'r' or band == 'g': # juntar
		ccdnum = 1
		while ccdnum <= 62:
			print '\t ccd %i/62'%ccdnum
			if ccdnum == 61: # el 61 no esta
				ccdnum = ccdnum+1
			else:
				# RR20_all_r_ccd10.fits
				archivo = 'RR%i/RR%i_all_%s_ccd%i.fits'%(field,field,band,ccdnum)
				table = Table.read(indir+'/'+archivo, format='fits')
				df = table.to_pandas()
	
				df = df[['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1


				#df.sort_values(by='MJDObs')
				#diff_epochs = df.MJDObs.unique()
	
				if ccdnum == 1:
					stacked_df = df
	
				if ccdnum > 1 :
					stacked_df = pd.concat([stacked_df, df], ignore_index=True)

				del df

				ccdnum = ccdnum+1

	stacked_df.to_csv(outdir+"/RR%i_all_%s.csv"%(field, band))

	del stacked_df

	field = field+1





 
