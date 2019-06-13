#!/usr/bin/python2.7

from __future__ import division
import matplotlib
matplotlib.use('Agg')
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
import P4J
from astropy.stats import bootstrap
from scipy.stats import  gumbel_r 
from astropy.stats import sigma_clip
import pandas as pd

from astroML.time_series import lomb_scargle



# ESTE PROGRAMA PLOTEA PERIODOGRAMAS PARA ALGUNAS ESTRELLAS, INCLUYENDO SUS 3 PERIODOS MAS SIGNIFICANTES. 




#indir = "/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Selected/Selected_si"
indir = "/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Selected/porseleccionar_si"
main_archivo = "01_with_3rdperiod"

outdir = indir+'/show3period'

mainfile = np.loadtxt( indir+'/'+main_archivo , dtype='str')
print mainfile
largo = len(mainfile)

i=0
while i<largo:

	archivo = mainfile[i]

	print i+1, '/', largo, '\t\t', archivo

	df = pd.read_csv(indir+'/'+archivo)

	mag = df['mag_rel_cal'].as_matrix()
	err = df['mag_rel_cal_err'].as_matrix()
	mjd = df['MJDObs'].as_matrix()
	
	my_per = P4J.periodogram(method='QMIEU')
	my_per.set_data(mjd, mag, err)
	my_per.frequency_grid_evaluation(fmin=0.5, fmax=10.0, fresolution=0.001)  # frequency sweep parameters
	my_per.finetune_best_frequencies(fresolution=0.0001, n_local_optima=4)
	freq, per = my_per.get_periodogram()
	fbest, pbest = my_per.get_best_frequencies() # Return best n_local_optima frequencies
	
	bestperiod = 1./fbest[0]
	bestperiod2 = 1./fbest[1]
	bestperiod3 = 1./fbest[2]
	bestperiod4 = 1./fbest[3]
	
	pbest_bootstrap = np.zeros(shape=(100, 2))
	for index in range(pbest_bootstrap.shape[0]):
		P = np.random.permutation(len(mjd))
		my_per.set_data(mjd, mag[P], err[P])
		my_per.frequency_grid_evaluation(fmin=0.0, fmax=4.0, fresolution=1e-3)
		my_per.finetune_best_frequencies(fresolution=1e-4, n_local_optima=pbest_bootstrap.shape[1])
		_, pbest_bootstrap[index, :] = my_per.get_best_frequencies()
	                                
		param = gumbel_r.fit(pbest_bootstrap.ravel())
		rv = gumbel_r(loc=param[0], scale=param[1])
		x = np.linspace(rv.ppf(0.001), rv.ppf(0.999), 100)
	                                
	p_vals = [0.01, 0.05, 0.08]
	sig1 = rv.ppf(1.-p_vals[0])
	sig5 = rv.ppf(1.-p_vals[1])
	sig8 = rv.ppf(1.-p_vals[2])
	
	bestpower = pbest[0]
	bestpower2 = pbest[1]
	bestpower3 = pbest[2]
	bestpower4 = pbest[3]
	
	
	
	print bestperiod, bestperiod2, bestperiod3, bestperiod4
	
	
	fig = plt.figure(figsize=(8, 8)) 
	# Bootstrap
	ax = fig.add_subplot(3, 2, 1)
	#_ = ax.hist(pbest_bootstrap.ravel(), density=True, bins=20, alpha=0.2, label='Peak\'s bootstrap')
	_ = ax.hist(pbest_bootstrap.ravel(), bins=20, alpha=0.2, label='Peak\'s bootstrap')
	ax.plot(x, rv.pdf(x), 'r-', lw=5, alpha=0.6, label='Gumbel PDF')
	ymin, ymax = ax.get_ylim()
	ax.plot([pbest[0], pbest[0]], [ymin, ymax], '-', linewidth=4, alpha=0.3, color='b', label="Per1 power")
	ax.plot([pbest[1], pbest[1]], [ymin, ymax], '-', linewidth=4, alpha=0.3, color='r', label="Per2 power")
	for p_val in p_vals:
		ax.plot([rv.ppf(1.-p_val), rv.ppf(1.-p_val)], [ymin, ymax], '--', linewidth=4, alpha=0.3, label=str(p_val))
	ax.set_ylim([ymin, ymax])
	plt.xlabel('Periodogram value')
	plt.title('Significance Test'); plt.legend(loc=1)
	
	# Periodogram
	ax = fig.add_subplot(3, 2, 2)
	ax.plot(freq, per, color='k')
	ymin, ymax = ax.get_ylim()
	ax.plot([fbest[0], fbest[0]], [ymin, ymax], linewidth=8, alpha=0.4)
	ax.plot([fbest[1], fbest[1]], [ymin, ymax], linewidth=8, alpha=0.4, color='red')
	ax.plot([fbest[2], fbest[2]], [ymin, ymax], linewidth=8, alpha=0.4, color='green')
	ax.plot([fbest[3], fbest[3]], [ymin, ymax], linewidth=8, alpha=0.4, color='grey')
	xmin, xmax = ax.get_xlim()
	for p_val in p_vals:
		ax.plot([xmin, xmax], [rv.ppf(1.-p_val), rv.ppf(1.-p_val)], '--', linewidth=4, alpha=0.5, label=str(p_val))
	ax.set_ylim([ymin, ymax])
	ax.set_xlabel('Frequency [1/MJD]')
	ax.set_ylabel('QMI Periodogram')
	plt.title('Periodogram'); plt.legend(title='p-val', loc=1);
	plt.grid()
	
	# best period
	ax = fig.add_subplot(3, 2, 3)
	phase = np.mod(mjd, 1.0/fbest[0])*fbest[0]
	idx = np.argsort(phase)
	ax.errorbar(np.concatenate([np.sort(phase)-1, np.sort(phase), np.sort(phase)+1.0]), 
		np.concatenate([mag[idx], mag[idx], mag[idx]]),
		np.concatenate([err[idx], err[idx], err[idx]]), fmt='.k', ecolor='red', markersize=12, elinewidth=1.5)
	plt.gca().invert_yaxis()
	plt.title('Best period     %0.5f [d]'%(1.0/fbest[0]))
	ax.set_xlabel('Phase @ %0.5f [1/d]' %fbest[0])
	ax.set_ylabel('Magnitude')
	plt.grid()
	plt.tight_layout()
		
	# 2nd best period
	ax = fig.add_subplot(3, 2, 4)
	phase2 = np.mod(mjd, 1.0/fbest[1])*fbest[1]
	idx2 = np.argsort(phase2)
	ax.errorbar(np.concatenate([np.sort(phase2)-1, np.sort(phase2), np.sort(phase2)+1.0]), 
		np.concatenate([mag[idx2], mag[idx2], mag[idx2]]),
		np.concatenate([err[idx2], err[idx2], err[idx2]]), fmt='.k', ecolor='red', markersize=12, elinewidth=1.5)
	plt.gca().invert_yaxis()
	plt.title('Best period2     %0.5f [d]'%(1.0/fbest[1]))
	ax.set_xlabel('Phase @ %0.5f [1/d]' %fbest[1])
	plt.grid()
	plt.tight_layout()
	
	# 3rd best period
	ax = fig.add_subplot(3, 2, 5)
	phase3 = np.mod(mjd, 1.0/fbest[2])*fbest[2]
	idx3 = np.argsort(phase3)
	ax.errorbar(np.concatenate([np.sort(phase3)-1, np.sort(phase3), np.sort(phase3)+1.0]), 
		np.concatenate([mag[idx3], mag[idx3], mag[idx3]]),
		np.concatenate([err[idx3], err[idx3], err[idx3]]), fmt='.k', ecolor='red', markersize=12, elinewidth=1.5)
	plt.gca().invert_yaxis()
	plt.title('Best period3     %0.5f [d]'%(1.0/fbest[2]))
	ax.set_xlabel('Phase @ %0.5f [1/d]' %fbest[2])
	plt.grid()
	plt.tight_layout()

	# 4th best period
	ax = fig.add_subplot(3, 2, 6)
	phase4 = np.mod(mjd, 1.0/fbest[3])*fbest[3]
	idx4 = np.argsort(phase4)
	ax.errorbar(np.concatenate([np.sort(phase4)-1, np.sort(phase4), np.sort(phase4)+1.0]), 
		np.concatenate([mag[idx4], mag[idx4], mag[idx4]]),
		np.concatenate([err[idx4], err[idx4], err[idx4]]), fmt='.k', ecolor='red', markersize=12, elinewidth=1.5)
	plt.gca().invert_yaxis()
	plt.title('Best period4     %0.5f [d]'%(1.0/fbest[3]))
	ax.set_xlabel('Phase @ %0.5f [1/d]' %fbest[3])
	plt.grid()
	plt.tight_layout()
	
	# Save
	fig.savefig( outdir+'/'+ archivo.replace('obj', 'starpngGLS').replace('.csv','.png'))
	plt.show(fig)
	plt.close(fig)
	
	
	i=i+1
	print '\n\n'
	





















