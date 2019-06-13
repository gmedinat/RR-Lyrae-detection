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
from pylab import *
from numpy import nan
import random
import subprocess
import lomb
import pyfits as fits
from astroML.time_series import lomb_scargle

#-juntador: lee todos los elementos en una carpeta (como 'Estrellas_Gen') y tambien los elementos de otra subcarpeta (como 'Candidatos'), que solo tenga imagenes y donde falten los .dat de esas imagenes. Luego los copia desde la primera carpeta a la segunda.

#dir1 = "Estrellas_Gen"
#dir2 = "Candidatos2.0"
#dir1 = "Estrellas_Gen"
#dir2 = "Candidatos"
dir1 = "/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Selected"
dir2 = "/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Selected/porseleccionar_si"


#files = glob.glob("/run/media/gmedina/TOSHIBA EXT/DECam2018/datos_leftraru/Estrellas_Gen_All_Periods_P4J/Selected/*.csv") 

dir1files = os.listdir(dir1)
dir2files = os.listdir(dir2)

contador = 0
for s2 in dir2files:
	#archivo   star_36.2134_2576.2747_26_S18.dat
	#imagen s2   starpngGLS_36.2134_2576.2747_26_S18.png
	print contador
	s1 = s2.replace('starpngGLS','obj').replace('png', 'csv')
	#os.system( 'cp -v %s/%s %s/'%( dir1, s1, dir2) ) # copy from   source  to  destiny
	os.system( 'mv %s/%s %s/'%( dir1, s1, dir2) ) # copy from   source  to  destiny
	
	contador = contador+1

print "FIN!"


