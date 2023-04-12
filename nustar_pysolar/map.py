import numpy as np
import sunpy.map
import sunpy.sun
import sunpy.coordinates

import astropy.time
import astropy.units as u


import logging


def make_sunpy(evtdata, hdr, exp_time=0, on_time=0, norm_map=False, 
               shevt_xy=(0,0)*u.arcsec, ssw_legacy=False):
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
    
    shevt_xy: 2 element array of x and y shift to apply to 
    	evtdata before making the map (does to nearest pixel),
    	defaults to [0,0], so no shift

    ssw_legacy: Using evt converted to S/C via older sswidl code,
            which might need to specify CUNIT1 and CUNIT2. Defaults to "False"
    
    Returns
    -------

    nustar_map:
        A sunpy map object
    
    """
    
    # Parse Header keywords
    for field in list(hdr.keys()):
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
                xval = field[5:8]
            if hdr[field] == 'Y':
                yval = field[5:8]

    min_x = hdr[f'TLMIN{xval}']
    min_y = hdr[f'TLMIN{yval}']
    max_x = hdr[f'TLMAX{xval}']
    max_y = hdr[f'TLMAX{yval}']

    xunit = u.Unit(hdr[f'TCUNI{xval}'])
    delx = abs(hdr[f'TCDLT{xval}']) * xunit

    yunit = u.Unit(hdr[f'TCUNI{yval}'])
    dely = abs(hdr[f'TCDLT{yval}']) * yunit

    x = evtdata['X'][:]
    y = evtdata['Y'][:]

    # Apply a shift to the data - default is 0
    x = x + round( (shevt_xy[0] / delx).decompose().value )
    y = y + round( (shevt_xy[1] / dely).decompose().value )

    met = evtdata['TIME'][:] * u.s
    mjdref = hdr['MJDREFI']

    mid_obs_time = astropy.time.Time(mjdref*u.d+met.mean(), format = 'mjd')
    sta_obs_time = astropy.time.Time(mjdref*u.d+met.min(), format = 'mjd')

    # Get the exposure and ontimes, just a numbers not units of seconds
    if (exp_time == 0):
        exp_time = hdr['EXPOSURE']
    if (on_time == 0):
        on_time = hdr['ONTIME']
        
    # Use the native binning for now
    # Assume X and Y are the same size
    resample = 1.0
    scale = delx * resample
    bins = (max_x - min_x) / (resample)

    # numpy error that histogram2d needs integer number of bins -> fixed  
    H, yedges, xedges = np.histogram2d(y, x, bins=int(bins), range=[[min_y,max_y], [min_x, max_x]])

    # Normalise the data with the exposure (or live) time?
    if norm_map is True:
        H = H / exp_time
        pixluname = 'DN/s'
    else:
        pixluname = 'DN'

    # Working with older sswidl converted evt?
    if ssw_legacy is True:
        cunit = "arcsec"
    else:
        cunit = scale.unit.to_string()

    dict_header = {
        # Change to the start time of the obs, not the mid_time?
        "DATE-OBS": sta_obs_time.iso,#mid_obs_time.iso,
        "EXPTIME": exp_time,
        "ONTIME": on_time,
        "CDELT1": scale.value,
        "NAXIS1": bins,
        "CRVAL1": 0.,
        "CRPIX1": bins*0.5,
        "CUNIT1": cunit,
        "CTYPE1": "HPLN-TAN",
        "CDELT2": scale.value,
        "NAXIS2": bins,
        "CRVAL2": 0.,
        "CRPIX2": bins*0.5 + 0.5,
        "CUNIT2": cunit,
        "CTYPE2": "HPLT-TAN",
        "PIXLUNIT": pixluname,
        "DETECTOR":"NuSTAR",
        "HGLT_OBS": sunpy.coordinates.sun.B0(mid_obs_time).value,
        "HGLN_OBS": 0,
        # previous form (i.e. -4d41m42.2166s, instead of -4.695060153501252) save to fits crashed
        "RSUN_OBS": sunpy.coordinates.sun.angular_radius(mid_obs_time).value,
        "RSUN_REF": sunpy.sun.constants.radius.value,   
        # again need in m, not AU or save to fits crashed
        "DSUN_OBS": sunpy.coordinates.sun.earth_distance(mid_obs_time).to_value('m')
    }

    header = sunpy.util.MetaDict(dict_header)
    nustar_map = sunpy.map.Map(H, header)
    
    return nustar_map

