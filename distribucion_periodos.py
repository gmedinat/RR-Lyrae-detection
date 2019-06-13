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
from astropy.io import fits 
import pandas as pd
from astropy.table import Table
import astropy
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy import stats

save = True
show = True

#indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J" 
#files = os.listdir("/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J")
indir  = "/run/media/gmedina/TOSHIBA EXT/borrar" 
files = os.listdir("/run/media/gmedina/TOSHIBA EXT/borrar")

i = 1
for f in files:
	#break
	if f[len(f)-4:]=='.csv': 
		
		print i, '/', len(files)/2.
		print '\t\t', f

		# obj_1_N1_306.1172403_-40.0546931_r.cvs
	
		#,RA_deg,DEC_deg,coord_ra,coord_dec,object,
		#Field,CCD,ccdnum,MJDObs,Obs,AIRMASS,SEEING,EXPTIME,
		#base_PsfFlux_flux,base_PsfFlux_fluxSigma,
		#a,a_err,b,b_err,k,k_err,mag_r,zp_rel,zp_rel_err,
		#mag_rel_cal,mag_rel_cal_err,period1,period2

		df = pd.read_csv(indir+'/'+f)

		RA_deg = df.iloc[0]['RA_deg']
		DEC_deg = df.iloc[0]['DEC_deg']
		coord_ra = df.iloc[0]['coord_ra']
		coord_dec = df.iloc[0]['coord_dec']

		Field = df.iloc[0]['Field']
		CCD = df.iloc[0]['CCD']
		ccdnum = df.iloc[0]['ccdnum']

		mean_mag = df['mag_rel_cal'].mean()

		period1 = df.iloc[0]['period1']
		period2 = df.iloc[0]['period2']


		dat = np.array([[RA_deg,DEC_deg,coord_ra,coord_dec,Field,CCD,ccdnum,mean_mag,period1,period2]])
		#if Matriz == None:
		if i == 1:
			Matriz = dat
		else:
			Matriz = np.vstack([Matriz, dat])
		del df
		i = i+1


#guardar la matriz en la carpeta en un archivo de texto, de modo que matchPY2magPrevGen.py la pueda leer
#np.savetxt( "/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J/001_period_distribution.dat", Matriz, fmt='%s',delimiter="   ", header='RA_deg \t DEC_deg \t coord_ra \t coord_dec \t Field \t CCD \t ccdnum \t mean_mag \t period1 \t period2')

np.savetxt( "./001_period_distribution_CHEQUEO.dat", Matriz, fmt='%s',delimiter="   ", header='RA_deg \t DEC_deg \t coord_ra \t coord_dec \t Field \t CCD \t ccdnum \t mean_mag \t period1 \t period2')

#Matriz = np.loadtxt('001_period_distribution_CHEQUEO.dat', dtype='str' )

P1 = Matriz[:,8]
P2 = Matriz[:,9]

P1 = map(float,P1)
P2 = map(float,P2)

print stats.mode(P1)
print stats.mode(P2)

fig, ax = plt.subplots()


#plt.hist(P1, bins=18)  # arguments are passed to np.histogram 
plt.hist(P1, bins=200)  # arguments are passed to np.histogram 

ax.set_xlim(0.1, 1)
#ax.set_ylim(0, 53)
#ax.set_xticks([16,17,18,19,20,21,22])
#ax.set_yticks([5,10,15,20,25,30,35,40,45,50])

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('Period 1', fontsize=20)
plt.ylabel('N', fontsize=20)

plt.tight_layout()
if save == True:
	#plt.savefig('histogram_P1.eps')
	#plt.savefig('histogram_P1.png')
	plt.savefig('histogram_P4J_P1.png')

if show == True:
	plt.show()








