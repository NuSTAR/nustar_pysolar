import astropy.time
import astropy.units as u
from astropy.coordinates import get_sun
import sunpy.map
from sunpy import sun

import numpy as np



def _get_sun_pos(met):
	""" Convenience wrapper module for getting the solar position
	and north pole angle at a given time.
	"""
	
	sun_time = astropy.time.Time(met, format = 'mjd')
	astro_sun_pos = get_sun(sun_time)

	# Get the center of the Sun, and assign it degrees.
	sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg

	# Solar NP roll angle:
	sun_np=sun.solar_north(t=sun_time).cgs
	return sun_pos, sun_np;


def _xy_to_radec(evtdata, hdr):
	# Parse the header information
	for field in hdr.keys():
		if field.find('TYPE') != -1:
			if hdr[field] == 'X':
				print(hdr[field][5:8])
				xval = field[5:8]
			if hdr[field] == 'Y':
				print(hdr[field][5:8])
				yval = field[5:8]

	ra_ref = hdr['TCRVL'+xval]*u.deg
	x0 = hdr['TCRPX'+xval]
	delx = hdr['TCDLT'+xval] * u.deg

	dec_ref = hdr['TCRVL'+yval]*u.deg
	y0 = hdr['TCRPX'+yval]
	dely = hdr['TCDLT'+yval]*u.deg


	# Make local copies for convenience
	x = evtdata['X']
	y = evtdata['Y']

	# Convert the NuSTAR epoch times to MJD
	mjdref=hdr['MJDREFI']
	met = evtdata['TIME']*u.s + mjdref*u.d

#	time = astropy.time.Time(mjdref*u.d+met, format = 'mjd')

	# Convert X and Y to RA/dec
	ra_x = ra_ref + (x - x0) * delx / np.cos(dec_ref)
	dec_y = dec_ref + (y - y0) * dely

	return ra_x, dec_y, met


def _delta_solar(ra_x, dec_y, met, **kwargs):

	# How often you want to update the solar ephemeris:
	tStep=kwargs.get('tStep', 5.0)
	tStep = tStep * u.s

	# How many events do you want to do?
	maxEvt=kwargs.get('maxEvt', len(ra_x))

	# Keep last time we updated things
	last_met = met[0] - tStep * 2.
	last_i = 0

	sun_x = np.zeros_like(ra_x)
	sun_y = np.zeros_like(dec_y)
	for i in np.arange(len(ra_x)):
		if( (met[i] - last_met) > tStep ):
			(sun_pos, sun_np) = _get_sun_pos(last_met)
			last_met = met[i]
			# Rotation matrix for a counter-clockwise rotation since we're going
			# back to celestial north from solar north
			rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
								  [-np.sin(sun_np),  np.cos(sun_np)]])

		# Diagnostics
		#         di = (i -last_i)
		#        print("Updating Sun position...")
		#         if di > 0:
		#             print(i, di)
		#             dt = toc()
		#             tic()
		#             last_i = i
		#             print("Time per event: ",dt / float(di) )
		# From here on we do things for every photon:

		ph_pos = np.array([ra_x[i].value, dec_y[i].value]) * u.deg
		offset = ph_pos - sun_pos
	   
		# Project the offset onto the Sun
		delta_offset = ((np.dot(offset, rotMatrix)).to(u.arcsec))

		# Account for East->West conversion for +X direction in heliophysics coords
		delta_offset = delta_offset*[-1., 1.]

		sun_x[i] = delta_offset[0]
		sun_y[i] = delta_offset[1]
		if (i>maxEvt):
			break
		
	return sun_x, sun_y


def to_solar(evtdata, hdr, **kwargs):

	# Convert events to RA/dec coordinates
	(ra_x, dec_y, met) = _xy_to_radec(evtdata, hdr)
	
	# Conver to solar coordinates
	(sun_x, sun_y) = _delta_solar(ra_x, dec_y, met, **kwargs)


	# Parse the header information to get the native bin size
	for field in hdr.keys():
		if field.find('TYPE') != -1:
			if hdr[field] == 'X':
				print(hdr[field][5:8])
				xval = field[5:8]
			if hdr[field] == 'Y':
				print(hdr[field][5:8])
				yval = field[5:8]
	delx = -1.0 * hdr['TCDLT'+xval] * u.deg
	dely = hdr['TCDLT'+yval]*u.deg


	# Make output variaibles:
	tbldata=evtdata.copy
	outhdr = hdr.copy

	#  change to 0-3000 pixels using the same pixel size.
	maxX = 3000
	maxY = 3000
	x0 = maxX / 2.
	y0 = maxY / 2.

	# Make image coordinates
	out_sun_x=(1.0)*(sun_x / delx) + x0
	out_sun_y=(sun_y / dely) + y0
	tbldata['X'] = out_sun_x
	tbldata['Y'] = out_sun_y

	# Update header information
	outhdr['TCRVL'+xval] = 0.
	outhdr['TCRPX'+xval] = x0
	outhdr['TCDLT'+xval] = 1.0 * delx.to(u.arcsec).value
	outhdr['TLMAX'+xval] = maxX
	outhdr['TCRVL'+yval] = 0.
	outhdr['TCRPX'+yval] = x0
	outhdr['TCDLT'+yval] = dely.to(u.arcsec).value
	outhdr['TLMAX'+yval] = maxY

	return tbldata, outhdr	
	
	

	
	








