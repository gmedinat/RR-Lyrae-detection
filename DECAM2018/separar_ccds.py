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

paridad = 'pares'
#paridad = 'impares'

if paridad == 'pares':
	indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/Fields_with_uncalibrated_magnitudes/pares" 
	outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/Fields_with_uncalibrated_magnitudes/pares/PerCCD" 

if paridad == 'impares':
	indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/Fields_with_uncalibrated_magnitudes/impares" 
	outdir = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/Fields_with_uncalibrated_magnitudes/impares/PerCCD" 

print '\n\n'

if paridad == 'pares':
	field = 24
if paridad == 'impares':
	field = 23

#falta el 23

while field <= 24 :

	print '\n Field ', field, '\n'	

	archivo = "RR%s_%s.csv"%(field, band)

	#file:///run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Fields_with_uncalibrated_magnitudes/RR1_g.csv
	
	df = pd.read_csv(indir+'/'+archivo)

	ccdnum_diferentes = df.ccdnum.unique()
	largo_ccdnum_diferentes = len(ccdnum_diferentes)


	#i = int(sys.argv[1])
	i = 0
	while i < largo_ccdnum_diferentes: 
		df_toFilter = df.loc[df['ccdnum'] == ccdnum_diferentes[i]]
		
		print '\t\t\t Saving ccdnum ', ccdnum_diferentes[i]
		df_toFilter.to_csv(outdir+"/RR%i_ccdnum%i_%s.csv"%(field, ccdnum_diferentes[i], band))

		del df_toFilter
		i = i+1

	del df

	field = field+2

	print '\n ------------------------------------------------------------------------------- \n'











 
