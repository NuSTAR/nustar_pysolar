import numpy as np
import sunpy.map
from sunpy import sun

import astropy.time
import astropy.units as u


import logging


def make_sunpy(evtdata, hdr):
	""" Make a sunpy map based on the NuSTAR data.
	
	Parameters
	----------
    evtdata: FITS data structure
		This should be an hdu.data structure from a NuSTAR FITS file.

	hdr: FITS header containing the astrometric information
	
	Returns
    -------

	nustar_map:
		A sunpy map objecct
	
	"""

	# Parse Header keywords
	for field in hdr.keys():
		if field.find('TYPE') != -1:
			if hdr[field] == 'X':
				print(hdr[field][5:8])
				xval = field[5:8]
			if hdr[field] == 'Y':
				print(hdr[field][5:8])
				yval = field[5:8]
		
	min_x= hdr['TLMIN'+xval]
	min_y= hdr['TLMIN'+yval]
	max_x= hdr['TLMAX'+xval]
	max_y= hdr['TLMAX'+yval]

	delx = hdr['TCDLT'+xval]
	dely = hdr['TCDLT'+yval]

	x = evtdata['X'][:]
	y = evtdata['Y'][:]
	met = evtdata['TIME'][:]*u.s
	mjdref=hdr['MJDREFI']
	mid_obs_time = astropy.time.Time(mjdref*u.d+met.mean(), format = 'mjd')

	# Use the native binning for now

	# Assume X and Y are the same size
	resample = 1.0
	scale = delx * resample
	bins = (max_x - min_x) / (resample)

	H, yedges, xedges = np.histogram2d(y, x, bins=bins, range = [[min_y,max_y], [min_x, max_x]])


	dict_header = {
	"DATE-OBS": mid_obs_time.iso,
	"CDELT1": scale,
	"NAXIS1": bins,
	"CRVAL1": 0.,
	"CRPIX1": bins*0.5,
	"CUNIT1": "arcsec",
	"CTYPE1": "HPLN-TAN",
	"CDELT2": scale,
	"NAXIS2": bins,
	"CRVAL2": 0.,
	"CRPIX2": bins*0.5 + 0.5,
	"CUNIT2": "arcsec",
	"CTYPE2": "HPLT-TAN",
	"HGLT_OBS": 0,
	"HGLN_OBS": 0,
	"RSUN_OBS": sun.solar_semidiameter_angular_size(mid_obs_time).value,
	"RSUN_REF": sun.constants.radius.value,
	"DSUN_OBS": sun.sunearth_distance(mid_obs_time).value
	}
	# For some reason the DSUN_OBS crashed the save...

	header = sunpy.map.MapMeta(dict_header)

	nustar_map = sunpy.map.Map(H, header)
	
	return nustar_map


