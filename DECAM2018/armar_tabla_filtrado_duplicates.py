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

files = glob.glob("/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Color_cut_test/*r.csv") 

largo = len(files)

Matriz = None
i=0
for f in files: 

	if i%500==0:
		print i, '/', largo

	print f

	df = pd.read_csv(f)

	RA = df.iloc[0]['RA_deg']
	DEC = df.iloc[0]['DEC_deg']
	Field = df.iloc[0]['Field']
	CCD = df.iloc[0]['CCD']
	CCDnum = df.iloc[0]['ccdnum']
	object_id = df.iloc[0]['object']
	N_r = df.iloc[0]['N_r'] #
	N_g = df.iloc[0]['N_g'] #
	mean_mag_rel_cal_r = df.iloc[0]['mean_mag_rel_cal_r']
	mean_mag_rel_cal_g = df.iloc[0]['mean_mag_rel_cal_g']
	g_r = mean_mag_rel_cal_g-mean_mag_rel_cal_r #
	g_r_distance = df.iloc[0]['g_r_distance'] #

	period1 = df.iloc[0]['period1']
	period2 = df.iloc[0]['period2']
	amplitud_total = abs( df['mag_rel_cal'].max() - df['mag_rel_cal'].min() )


	MEAN_MAG_REL_CAL_ERR_r = df['mag_rel_cal_err'].mean()

	dat = np.array([[RA, DEC, Field, CCD, CCDnum, object_id, period1, period2, N_r, N_g, amplitud_total, mean_mag_rel_cal_r, mean_mag_rel_cal_g, MEAN_MAG_REL_CAL_ERR_r, g_r, g_r_distance]])

	if i==0:
		Matriz = dat
	else:
		Matriz=np.vstack([Matriz, dat])



	i=i+1

np.savetxt( "/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Color_cut_test/001Tabla.dat", Matriz, fmt='%s',delimiter="   ", header='RA \t DEC \t Field \t CCD \t CCDnum \t object_id \t period1 \t period2 \t N_r \t N_g \t amplitud_total \t mean_mag_rel_cal_r \t mean_mag_rel_cal_g \t MEAN_MAG_REL_CAL_ERR_r \t mean_mag_rel_cal_g-mean_mag_rel_cal_r \t g_r_distance')
np.savetxt( "/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Color_cut_test/001Tabla.csv", Matriz, fmt='%s',delimiter=",")
#RA,DEC,Field,CCD,CCDnum,object_id,period1,period2,N_r,N_g,amplitud_total,mean_mag_rel_cal_r,mean_mag_rel_cal_g,MEAN_MAG_REL_CAL_ERR_r,mean_mag_rel_cal_g-mean_mag_rel_cal_r,g_r_distance









