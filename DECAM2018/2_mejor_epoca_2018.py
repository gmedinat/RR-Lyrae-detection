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
import timeit

start = timeit.default_timer()

#filtro = 'g'
filtro = 'r'

campo = 1

Matriz = None

while campo <=24:

	indir  = "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/RR%s"%campo
	files = os.listdir("/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/RR%s"%campo)
	
	print '\nField %s *********************************'%campo

	i = 1
	for f in files:

		#RR23_18_58231.42175202_r.csv
		#RR23_1_58229.22433551_r.csv
		#RR2_18_58229.22433551_r.csv
		#RR2_1_58229.22433551_r.csv
		#RR23_all_g.csv
		#RR23_all_g.fits
		#RR23_all_r.csv
		
		#print f, f[len(f)-9:len(f)-6]
		if f[len(f)-5:]=='%s.csv'%filtro and f[len(f)-9:len(f)-6] != 'all':
			
			if f[3]=='_' and f[5]=='_': #RR2_1_58229.22433551_r.csv
				epoca = f[4]
				#busqueda_guion
				j = 6
				while f[j]!='_':
					j=j+1
				MJD = f[6:j]
			if f[4]=='_' and f[6]=='_': #RR23_1_58229.22433551_r.csv
				epoca = f[5]
				#busqueda_guion
				j = 7
				while f[j]!='_':
					j=j+1
				MJD = f[7:j]
			if f[3]=='_' and f[6]=='_': #RR2_18_58229.22433551_r.csv
				epoca = f[4:6]
				#busqueda_guion
				j = 7
				while f[j]!='_':
					j=j+1
				MJD = f[7:j]
			if f[4]=='_' and f[7]=='_': #RR23_18_58229.22433551_r.csv
				epoca = f[5:7]
				#busqueda_guion
				j = 8
				while f[j]!='_':
					j=j+1
				MJD = f[8:j]

			print '\t \t', f, '\t', epoca, MJD 
			#table = Table.read(indir+'/'+f, format='fits')
			#df = table.to_pandas()
			df = pd.read_csv(indir+'/'+f)

			xx=df['base_SdssShape_psf_xx']
			yy=df['base_SdssShape_psf_yy']
			df['rad'] = 0.263*np.sqrt(xx**2 + yy**2)   # where 0.263''/pix is the plate scale of DECam

			#MJD = df.loc[0, 'MJDObs']
			# sacar MJD del nombre del archivo mejor, para evitar problemas con el formato decimal
			
			

			mean = df['rad'].mean()
			
			dat = np.array([[MJD, epoca, mean]])
			#if Matriz == None:
			if i == 1:
				Matriz = dat
			else:
				Matriz = np.vstack([Matriz, dat])

			del df
			i = i+1

	#ordenar Matriz segun la 1ra columna
	Matriz=Matriz[Matriz[:,0].argsort()]

	#guardar la matriz en la carpeta en un archivo de texto, de modo que matchPY2magPrevGen.py la pueda leer
	np.savetxt( "/run/media/gmedina/TOSHIBA EXT/DECam2018/Data/mean_rad_Field_%s_%s.dat"%(campo, filtro), Matriz, fmt='%s',delimiter="   ", header='MJD \t Obs \t mean_rad')
			
	MJDs = Matriz[:,0].astype(np.float)
	means = Matriz[:,2].astype(np.float)

	argmini = np.argmin(means)

	fig, ax = plt.subplots()
	#ax.errorbar(Matriz[:,0], Matriz[:,2], fmt='ko', markersize=15, label='Field: %s'%campo)
	ax.errorbar(MJDs, means, fmt='ko', markersize=15, label='Field: %s'%campo)
	ax.errorbar(MJDs[argmini], means[argmini], fmt='ro', markersize=15)
	#ax.errorbar(u_g_hits, g_r_hits, xerr=u_g_hits_err, yerr=g_r_hits_err, fmt='ro',label='HiTS')
	#plt.legend(loc = 3, fontsize = 10, framealpha = 0.5)
	#ax.set_xlim(-0.5, 6.5)
	#ax.set_ylim(-0.45, 1.89)
	ax.grid(True)

	ax.set_xlabel('MJD', fontsize=20)
	ax.set_ylabel('rad', fontsize=20)
	ax.set_title('Field %s , %s-band      Min: %.4f (epoch %s)'%(campo, filtro, MJDs[argmini], Matriz[argmini,1]), fontsize=15)
	plt.xticks(fontsize=8)
	plt.yticks(fontsize=10)


#plt.title(r'{\fontsize{30pt}{3em}\selectfont{}{Mean WRFv3.5 LHF\r}{\fontsize{18pt}{3em}\selectfont{}(September 16 - October 30, 2012)}')


	ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

	plt.tight_layout()
	plt.savefig(indir+'/radius_plot_Field_%s_%s.png'%(campo,filtro))

	#plt.show()
	plt.close()
	#plt.savefig('color_color.eps')

	campo = campo+1 


stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))






