import numpy as np
import sunpy.map
import sunpy.sun
import sunpy.coordinates

import astropy.time
import astropy.units as u


import logging


def make_sunpy(evtdata, hdr, exp_time=0,on_time=0,norm_map=False):
    """ Make a sunpy map based on the NuSTAR data.
    
    Parameters
    ----------
    evtdata: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    hdr: FITS header containing the astrometric information

    Optional keywords
    
    exp_time: The exposure time (i.e. livetime, not on-time) no units. 
            If not given then taken from hdr
    on_time: The on-time (i.e. dwell) no units. 
            If not given then taken from hdr        

    norm_map: Normalise the map data by the exposure time (i.e. livetime), 
        giving map in units of DN/s. Defaults to "False" and units of DN
    
    Returns
    -------

    nustar_map:
        A sunpy map objecct
    
    """
    # Parse Header keywords
    for field in hdr.keys():
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
#                 print(hdr[field][5:8])
                xval = field[5:8]
            if hdr[field] == 'Y':
#                 print(hdr[field][5:8])
                yval = field[5:8]

    min_x= hdr['TLMIN'+xval]
    min_y= hdr['TLMIN'+yval]
    max_x= hdr['TLMAX'+xval]
    max_y= hdr['TLMAX'+yval]

    delx = abs(hdr['TCDLT'+xval])
    dely = abs(hdr['TCDLT'+yval])

    x = evtdata['X'][:]
    y = evtdata['Y'][:]
    met = evtdata['TIME'][:]*u.s
    mjdref=hdr['MJDREFI']

    mid_obs_time = astropy.time.Time(mjdref*u.d+met.mean(), format = 'mjd')
    sta_obs_time = astropy.time.Time(mjdref*u.d+met.min(), format = 'mjd')

    # Get the exposure and ontimes, just a numbers not units of seconds
    if (exp_time ==0):
        exp_time=hdr['EXPOSURE']
    if (on_time ==0):
        on_time=hdr['ONTIME']    

    # Use the native binning for now

    # Assume X and Y are the same size
    resample = 1.0
    scale = delx * resample
    bins = (max_x - min_x) / (resample)

#   numpy error that histogram2d needs integer number of bins -> fixed  
    H, yedges, xedges = np.histogram2d(y, x, bins=int(bins), range = [[min_y,max_y], [min_x, max_x]])

    #Normalise the data with the exposure (or live) time?
    if norm_map is True:
    	H=H/exp_time
    	pixluname='DN/s'
    else:
    	pixluname='DN' 

    dict_header = {
#    Change to the start time of the obs, not the mid_time?
    "DATE-OBS": sta_obs_time.iso,#mid_obs_time.iso,
    "EXPTIME": exp_time,
    "ONTIME": on_time,
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
    "PIXLUNIT": pixluname,
    "DETECTOR":"NuSTAR",
#    SunPy v1+ version
    "HGLT_OBS": sunpy.coordinates.sun.B0(mid_obs_time),
#    SunPy <v1 version        
#   "HGLT_OBS": sunpy.coordinates.get_sun_B0(mid_obs_time),
    "HGLN_OBS": 0,
#    SunPy v1+ version
    "RSUN_OBS": sunpy.coordinates.sun.angular_radius(mid_obs_time).value,
#    SunPy <v1 version
#    "RSUN_OBS": sun.solar_semidiameter_angular_size(mid_obs_time).value,
##        
#    sun.constants will stay in SunPy v1+ ? just sun.sun being depreciated?
#      https://docs.sunpy.org/en/stable/code_ref/sun.html#module-sunpy.sun.models
    "RSUN_REF": sunpy.sun.constants.radius.value,
###       
    # Assumes dsun_obs in m if don't specify the units, so give units
#    SunPy v1+ version        
    "DSUN_OBS": sunpy.coordinates.sun.earth_distance(mid_obs_time).value*u.astrophys.au
    #    SunPy <v1 version    
#     "DSUN_OBS": sunpy.coordinates.get_sunearth_distance(mid_obs_time).value*u.astrophys.au
    }
    # For some reason the DSUN_OBS crashed the save...
    
#    header = sunpy.map.MapMeta(dict_header)
    header = sunpy.util.MetaDict(dict_header)
    nustar_map = sunpy.map.Map(H, header)
    
    return nustar_map


