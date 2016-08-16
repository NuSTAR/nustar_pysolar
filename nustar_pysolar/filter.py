import numpy as np
import logging



def bad_pix(evtdata, fpm='A'):
	"""Do some basic filtering on known bad pixels.
	
	Parameters
	----------
    evtdata: FITS data structure
		This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"FPMA" | "FPMB"}
		Which FPM you're filtering on. Assumes A if not set.


	Returns
    -------

	goodinds: iterable
		Index of evtdata that passes the filtering.
	"""
    
	logging.basicConfig(filename='nustar_pysolar.log', level=logging.INFO)
	logging.info('Started')

	
	# Hot pixel filters
	
	# FPMA or FPMB
	
	if fpm.find('B') == -1 :
		logging.info("Filtering bad pixels for FPMA.")
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
		logging.info("Filtering bad pixels for FPMB.")
		pix_filter = np.invert( ( (evtdata['DET_ID'] == 0) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 24)) )


	inds = (pix_filter).nonzero()
	goodinds=inds[0]
	
	return goodinds
	
def by_energy(evtdata, energy_low=2.5, energy_high=10.):
	""" Apply energy filtering to the data.
	
	Parameters
	----------
	evtdata: FITS dat structure
		This should be an hdu.data structure from a NuSTAR FITS file.
		
	energy_low: float
		Low-side energy bound for the map you want to produce (in keV).
		Defaults to 2.5 keV.

    energy_high: float
		High-side energy bound for the map you want to produce (in keV).
		Defaults to 10 keV.
	"""		
	pilow = (energy_low - 1.6) / 0.04
	pihigh = (energy_high - 1.6) / 0.04
	pi_filter = ( ( evtdata['PI']>pilow ) &  ( evtdata['PI']<pihigh))
	inds = (pi_filter).nonzero()
	goodinds=inds[0]
	
	return goodinds
	
def gradezero(evtdata):
	""" Only accept counts with GRADE==0.
		
	Parameters
	----------
	evtdata: FITS dat structure
		This should be an hdu.data structure from a NuSTAR FITS file.
		
	Returns
    -------

	goodinds: iterable
		Index of evtdata that passes the filtering.
	"""

	# Grade filter
	
	grade_filter = ( evtdata['GRADE'] == 0)
	inds = (rade_filter).nonzerio()
	goodinds = inds[0]
	
	return goodinds


