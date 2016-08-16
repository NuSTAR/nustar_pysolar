import numpy as np
import logging



def filter_events(evtdata, fpm='A',energy_low=2.5, energy_high=10.):
	"""Do some basic filtering on known bad pixels.
	
	Parameters
	----------
    evtdata: FITS data structure
		This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"FPMA" | "FPMB"}
		Which FPM you're filtering on. Assumes A if not set.

    energy_low: float
		Low-side energy bound for the map you want to produce (in keV).

    energy_high: float
		High-side energy bound for the map you want to produce (in keV).

	Returns
    -------

	goodinds: iterable
		Index of evtdata that passes the filtering.
	"""
    
	logging.basicConfig(filename='nustar_pysolar.log', level=logging.INFO)
	logging.info('Started')

	pilow = (energy_low - 1.6) / 0.04
	pihigh = (energy_high - 1.6) / 0.04

	# Grade filter
	grade_filter = ( evtdata['GRADE'] == 0)
	pi_filter = ( ( evtdata['PI']>pilow ) &  ( evtdata['PI']<pihigh))
	
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


	inds = (grade_filter & pi_filter & pix_filter).nonzero()
	goodinds=inds[0]

	logging.info("Found: ", len(goodinds), " good counts.")

	logging.info('Finished')

	return goodinds


