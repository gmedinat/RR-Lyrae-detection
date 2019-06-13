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
import scipy
from pylab import *
from numpy import nan
import random
import subprocess
import pandas as pd
import glob

#outdir = '../'

result = None

filtro = 'r'
catalogo_referencia = 8

i=1
while i<=16:

	archivo_list = glob.glob("/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_%s/Field%s_%s_*_%s.csv"%(i,i,catalogo_referencia,filtro))
	#file:///run/media/gmedina/TOSHIBA EXT/DECam2017/Data/Field_1/Field1_8_57993.18054964_r.csv
	#if i < 10:
	#	nombre = 'newRADEC_catalog02_BlindA-P_0%i'%(i)
		#nombre = 'sdss%i_v4_GMedina.dat'%i
	#else: 
	#	nombre = 'newRADEC_catalog02_BlindA-P_%i'%(i)
		#nombre = 'sdss%i_v4_GMedina.dat'%i

	#print i, nombre

	#Field1_8_57993.18054964_r.csv

	archivo = archivo_list[0]
	print archivo

	df = pd.read_csv(archivo)

	if i == 1:
		big_df = df
	else:
		big_df = pd.concat([big_df,df])

	del df

	#data = np.loadtxt(indir+'/'+nombre, dtype='str')

	#if len(data)>0:
	#	data = np.reshape( data, (len(data), len(data[0,:])) )

	#	print len(data)

	#	if i == 1:

	#		result = data
	

	#	else:
			#result = pd.concat([result,data], axis = 1) #axis=1, join='inner'
			#result = np.reshape( result, (len(result), 1) )
	
			#result = result.append(data)
	#		result = np.vstack((result, data))
		
	#		print len(result)




	i = i+1

big_df.to_csv( '/run/media/gmedina/TOSHIBA EXT/DECam2017/Data/All_Fields_ref_%s.csv'%filtro )
del big_df

#np.savetxt( 'hipertabla_HiTS_newRADEC', result, fmt='%s %s %s %s %s %s',delimiter="   ", header='ra   dec   mag   emag   CCD   CCDnum')









