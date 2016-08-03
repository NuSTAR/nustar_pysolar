#!/usr/bin/env python

import sys, getopt, os
from os.path import *

from astropy.io import fits
import astropy.time
import astropy.units as u

from sunpy import sun
import sunpy.map

import numpy as np

scriptname='nustar_make_sunpy_map.py'



def main(argv):
	infile = ''	

	try:
		opts, args = getopt.getopt(argv,"hi:",["ifile="])

	except getopt.GetoptError:
		print(scriptname+'-i <inputfile>')
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print(scriptname+'-i <inputfile> -t <timestep in seconds>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			infile = arg


	print(scriptname)
	print("Processing: ", infile)


### Get the data from the FITS file.
	hdulist = fits.open(infile)

	evtdata=hdulist[1].data
	hdr = hdulist[1].header

	print("Loaded: ", len(evtdata['X']), " counts.")
	print("Effective exposure: ", hdr['EXPOSURE'], ' seconds')
	hdulist.close()
	

	# Do some basic filtering
	
	energy_low = 2.5
	energy_high = 10.
	pilow = (energy_low - 1.6) / 0.04
	pihigh = (energy_high - 1.60) /  0.04

	# Grade filter
	grade_filter = ( evtdata['GRADE'] == 0)
	pi_filter = ( ( evtdata['PI']>pilow ) &  ( evtdata['PI']<pihigh))

	# Hot pixel filters:
	# FPMA or FPMB
	if infile.find('B06') == -1 :
		print("Filtering bad pixels for FPMA.")
		pix_filter = np.invert( ( (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 16) & (evtdata['RAWY'] == 5) |
								(evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 22) |
								(evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 27) & (evtdata['RAWY'] == 6) |
								(evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 27) & (evtdata['RAWY'] == 21) |
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 22) & (evtdata['RAWY'] == 1) |
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 15) & (evtdata['RAWY'] == 3) |
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 5) & (evtdata['RAWY'] == 5) | 
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 22) & (evtdata['RAWY'] == 7) | 
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 16) & (evtdata['RAWY'] == 11) | 
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 18) & (evtdata['RAWY'] == 3) | 
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 4) | 
								(evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 25) & (evtdata['RAWY'] == 5) ) )
	else:
		print("Filtering bad pixels for FPMB.")
		pix_filter = np.invert( ( (evtdata['DET_ID'] == 0) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 24)) )


	inds = (grade_filter & pi_filter & pix_filter).nonzero()
	goodinds=inds[0]

	print("Found: ", len(goodinds), " good counts.")


	
	# Make the sunpy map and save it to the output file.
	
	# Header keywords
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

	x = evtdata['X'][goodinds]
	y = evtdata['Y'][goodinds]
	met = evtdata['TIME'][goodinds]*u.s

	# Conver the NuSTAR epoch times to astropy datetime objects
	mjdref=hdr['MJDREFI']
	mid_obs_time = astropy.time.Time(mjdref*u.d+met.mean(), format = 'mjd')

    
	# Use the native binning for now

	# Assume X and Y are the same size
	resample = 5.0
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
        "RSUN_REF": sun.constants.radius.value
    }

#        "DSUN_OBS": sun.sunearth_distance(mid_obs_time) * sun.constants.au.value
    
	header = sunpy.map.MapMeta(dict_header)
	nustar_map = sunpy.map.Map(H, header)

	# # Make the new filename:
	(sfile, ext)=splitext(infile)
	outfile=sfile+'_map.fits'
	print("Dumping to: ", outfile)

	# Remove output file if necessary
	if isfile(outfile):
		print(outfile, 'exists! Removing old version...')
		os.remove(outfile)

	nustar_map.save(outfile, filetype='auto')




if __name__ == "__main__":
   main(sys.argv[1:])
