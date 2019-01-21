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
#import lomb
import pyfits as fits
#from astroML.time_series import lomb_scargle
import glob
from astropy.stats import sigma_clip
import pandas as pd
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import timeit
start = timeit.default_timer()

obs_night_mode = 'mjd'
obs_night_mode = 'ut'

obs_night = 58497.0

#mode = 'show_airmass'
mode = 'show_phaseOnly'
#mode = 'show Low-MAX'

#times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']
#times = ['2019-01-14T00:00:00.00']
#t = Time(times, format='isot', scale='utc')
#print t.mjd  


#location          example from: http://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html
bear_mountain = EarthLocation(lat=41.3*u.deg, lon=-74*u.deg, height=390*u.m)
location = EarthLocation.of_site(u'Las Campanas Observatory')

#m33 = SkyCoord.from_name('M33') # ra dec in deg
#m33 = SkyCoord(152.481579691, 1.5366530542, unit='deg', frame='icrs') # frame: icrs, galactic, FK4, ...


#utcoffset = -4*u.hour  # Eastern Daylight Time
#time = Time('2012-7-12T23:00:00') - utcoffset
#m33altaz = m33.transform_to(AltAz(obstime=time,location=bear_mountain))
#m33altaz = m33.transform_to(AltAz(obstime=Time('2019-01-14T00:00:00'),location=bear_mountain))

#print m33

#midnight = Time('2019-01-14T00:00:00') - utcoffset
#midnight = Time('2019-01-14T00:00:00') - utcoffset
#delta_midnight = np.linspace(-2, 10, 12)*u.hour
#print delta_midnight
#frame_July13night = AltAz(obstime=midnight+delta_midnight, location=bear_mountain)
#m33altazs_July13night = m33.transform_to(frame_July13night)
#m33airmasss_July13night = m33altazs_July13night.secz
#frame_Jan14night = AltAz(obstime=midnight+delta_midnight, location=location)
#m33altazs_Jan14night = m33.transform_to(frame_Jan14night)
#m33airmasss_Jan14night = m33altazs_Jan14night.secz

#print m33altazs_Jan14night
#print 'airmasses: ', m33airmasss_Jan14night

#print EarthLocation.get_site_names()
#u'Cerro Pachon', u'Cerro Paranal', u'Cerro Tololo', u'Cerro Tololo Interamerican Observatory'
#u'La Silla Observatory', u'Large Binocular Telescope', u'Las Campanas Observatory'

print '\n\n'


archivo = 'results_fit_test_targets.dat'
# read table
data = np.loadtxt( './%s'%archivo , dtype='str')

ids = data[:,0]
periods = data[:,4]
mean_mags = data[:,5]
max_lights = data[:,9]

RAs  = data[:,1]
DECs = data[:,2]

largo = len(ids)

j = 0
while j < len(data):
	
	print j+1, '\t', ids[j]
	j=j+1

print '\n'

date = '2019-01-15'
midnight = Time('%sT00:00:00'%date)
nightend = Time('%sT09:30:20'%date)
delta_time = 30*u.min

#print dthr, dtmin, dtsec
#print dt/u.hr

#print dtsec.value, dtsec.value/60, dtsec.value%60
#print dtsec.value/3600
#dtt = dt - int(dthr)

print '(yyyy-mm-dd)\tUT\t\tMJD\t\t', 
j = 0
while j < len(data):
	
	print j+1, '\t',
	j=j+1
print '\n--------------------------------------------------------------------------------------------------------------------------'

print '           \t  \t\t   \t\t', 
j = 0
while j < len(data):
	
	print mean_mags[j], '\t',
	j=j+1
print '\n##############################################################################################################################'

sample_time = midnight
while sample_time < nightend:

	t = Time(sample_time, format='isot', scale='utc')

	dt = sample_time - midnight
	dthr = dt.to(u.hr)	
	dthrval = dthr.value
	hr = int(dthrval)
	m = int( (dthrval - hr)*60. )
	s = ( (dthrval - hr)*60. - m )*60.
	
	if 60 - s < 0.001:
		s = 0
		m = m+1

	if 60 - m < 0.001:
		m = 0
		hr = hr+1

	if s < 0.001:
		s = 0
	if m < 0.001:
		m = 0
	#dts = dt.to(u.s)
	#dtsval = dts.value
	#hr = dtsval/3600.
	#hr_str, m_str, s_str = hr, m, s

	if hr >= 10:
		hr_str = str(hr)
	else:
		hr_str = '0'+str(hr)
	if m >= 10:
		m_str = str(m)
	else:
		m_str = '0'+str(m)	
	if s >= 10:
		s_str = str(s)
	else:
		s_str = '0'+str(s)
	
	#print date, t, '%.2f'%(t.mjd), '\t\t', 
	print date, '\t', hr_str, m_str, s_str, '\t','%.2f'%(t.mjd), '\t', 

	i = 0
	while i<largo:

		RA  = RAs[i]
		DEC = DECs[i]

		obj_coords = SkyCoord(RA,DEC, unit='deg', frame='icrs') # frame: icrs, galactic, FK4, ...

		frame_AltAz = AltAz(obstime=sample_time, location=location)
		obj_AltAz = obj_coords.transform_to(frame_AltAz)
		obj_Airmass = obj_AltAz.secz


		ID = ids[i]
		epoch_max_light = float(max_lights[i])
		period = float(periods[i])

		diff = t.mjd - epoch_max_light
		if i != (largo-1):
			#print '%.2f'%(diff%period) , '\t',
			if mode == 'show_airmass':
				print '%.2f'%((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period )) , '_','%.1f'%obj_Airmass,'\t',

			if mode == 'show_phaseOnly':
				print '%.2f'%((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period )) , '\t',

			if mode == 'show Low-MAX':
				if obj_Airmass >= 1.4 or obj_Airmass < 1.:
					print 'Low' , '\t',
				elif ((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period ) < 0.10) or ((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period ) >= 0.90):
					print 'MAX' , '\t',
				else:
					print '%.2f'%((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period )) , '\t',

		else:
			#print '%.2f'%(diff%period) , '\t'
			if mode == 'show_airmass':
				print '%.2f'%((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period )) , '_','%.1f'%obj_Airmass,'\t'

			if mode == 'show_phaseOnly':
				print '%.2f'%((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period )) , '\t'

			if mode == 'show Low-MAX':
				if obj_Airmass > 1.3 or obj_Airmass < 1.:
					print 'Low' , '\t'
				elif ((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period ) < 0.10) or ((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period ) >= 0.90):
					print 'MAX' , '\t'
				else:
					print '%.2f'%((t.mjd-epoch_max_light)/period - int( (t.mjd-epoch_max_light)/period )) , '\t'


		#sys.exit(0)

		i = i+1

	#print ''

	#print sample_time
	sample_time = sample_time + delta_time


#print '\n'

#i = 0
#while i<largo:
#	ID = ids[i]
#	epoch_max_light = float(max_lights[i])
#	period = float(periods[i])

#	date = '2019-01-14'

#	print 'max light epoch: ', epoch_max_light , '\t period (d):', period , '\t to be observed at: ', obs_night
#	print 'julian date difference: ', obs_night - epoch_max_light

#	print '\n\n'
#	hr = 1
#	while hr <= 9:
#		if hr <=9:
#			str_hr = '0'+str(hr)
#		if hr > 9:
#			str_hr = str(hr)
#
#		times = ['%sT%s:00:00.00'%(date,str_hr)]
#		t = Time(times, format='isot', scale='utc')
#
#		diff = t.mjd - epoch_max_light
		#print t.mjd , '\t\t diff: ', diff, '\t modulo periodo: ', diff%period 
#		print date, str_hr, t.mjd, diff%period 
#		hr = hr+1

	


print '\n'





stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))
