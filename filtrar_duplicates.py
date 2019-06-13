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


f = "/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Color_cut_test/001Tabla.csv"
df = pd.read_csv(f)

largo = len(df)
print largo


index_concat = 1
i = 0
while True:
	if i%100==0:
		print '\t\t', i, '/', largo
	RA = df.iloc[i]['RA']
	DEC = df.iloc[i]['DEC']
	Field = df.iloc[i]['Field']
	object_id = df.iloc[i]['object_id']

	df['distancia_deg'] = (((df['RA']-float(RA))**2.+(df['DEC']-float(DEC))**2.)**0.5)

	df_in = df[(df['distancia_deg'] < 0.5/3600.)]

	#print df['distancia_deg']
	#print len(df_in), RA,DEC, object_id
	#print df_in
	#break

	if len(df_in) > 1:
		#df_sorted = df_in.sort(['N_r'])
		df_in_sorted = df_in.sort_values(by='N_r', ascending=0) 
		df_max = df_in_sorted[df_in_sorted['N_r'] == df_in_sorted.iloc[0]['N_r']] # contiene los que tengan mas N_r

		#print df_in_sorted
		#print df_max

		#print df.columns.values
		if len(df_max) > 1: # hay mas de 1 obj_id con el maximo de observaciones
			#asdsa #conservar solo el con menor error
		
			index_min = df_max['MEAN_MAG_REL_CAL_ERR_r'].idxmin()
			selected_df = df_max.loc[index_min]
			selected_df = selected_df.to_frame()
			#pd.DataFrame(df1).T
			selected_df = pd.DataFrame(selected_df)

			#selected_df = df_max.iloc[0]
			#df_max_sorted = df_max.sort_values(by='MEAN_MAG_REL_CAL_ERR_r', ascending=1) 
			#selected_df = df_in_sorted[df_in_sorted['N_r'] == df_in_sorted.iloc[0]['N_r']]

			#print selected_df.T
			#print df_max.shape, selected_df.shape

			selected_df = selected_df.T

			if index_concat == 1:
				big_df = selected_df
				index_concat = index_concat+1
			else:
				big_df = pd.concat([big_df,selected_df])
		else: 
			selected_df = df_max
			#selected_df = selected_df.T
			if index_concat == 1:
				big_df = selected_df
				index_concat = index_concat+1
			else:
				big_df = pd.concat([big_df,selected_df])
	else:
		selected_df = df_in
		#print selected_df
		if index_concat == 1:
			big_df = selected_df
			index_concat = index_concat+1
		else:
			big_df = pd.concat([big_df,selected_df])

	
	i = i+1

	if i == largo or i == -100:
		break

#filtro duplicates
#big_df_duplicates_removed = big_df
big_df_duplicates_removed = big_df.drop_duplicates(subset=['RA','DEC','object_id'], keep='first')

big_df_duplicates_removed.to_csv( "/run/media/gmedina/TOSHIBA EXT/DECam2017/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Color_cut_test/001Tabla_filtrada.csv" )





 
