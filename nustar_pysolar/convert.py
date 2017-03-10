import astropy.time
import astropy.units as u
from astropy.coordinates import get_sun
import sunpy.map
from sunpy import sun
import numpy as np


def _get_sun_pos(met):
    """
    Convenience wrapper module for getting the solar position
    and north pole angle at a given time.

    Synax:
    ----------
    sun_pos, sun_np = _get_sun_pos(met)

    Parameters
    ----------
    met: Time for the observation, given in MJD.

    Returns
    -------

    sun_pos:
        Index of evtdata that passes the filtering.
    """ 
    
    sun_time = astropy.time.Time(met, format = 'mjd')
    astro_sun_pos = get_sun(sun_time)

    # Get the center of the Sun, and assign it degrees.
    sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg

    # Solar NP roll angle:
    sun_np=sun.solar_north(t=sun_time).cgs
    return sun_pos, sun_np;

def _xy_to_radec(evtdata, hdr):
    """ Conversion function to go from X/Y coordinates
        in the FITS file to RA/Dec coordinates.
    """
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

#   time = astropy.time.Time(mjdref*u.d+met, format = 'mjd')

    # Convert X and Y to RA/dec
    ra_x = ra_ref + (x - x0) * delx / np.cos(dec_ref)
    dec_y = dec_ref + (y - y0) * dely

    return ra_x, dec_y, met


def _delta_solar(ra_x, dec_y, met, **kwargs):
    """ Function to compute the offsets from the center of
        the Sun as a function of time.

        Use the tStep argument to define how often you want
        to update the solar ephemeris. Default is every 5 seconds.

        Inputs: ra_x, dec_y, and met are all arrays that contain the
        RA, Dec, and time of arrival of each count, respectively.

        Outputs: sun_x, sun_y are the x and y values (in arcseconds)
        from the center of the Sun.

    """
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
    """ Main script to convert the events to solar coordinates.

        Inputs:
        --------
        evtdata: input event data from a NuSTAR Level2 FITS file.
        (i.e. nu*A06_cl.evt)
        
        hdr: The header from the same NuSTAR FITS file.
        
        Outputs:
        --------
        tbldata: Binary table with the positions converted from RA/Dec
            to heliocentric coordinates.
            
        outhdr: A new FITS header with the appropraite keywords adjusted for
            heliocentric coordinates.
            
        Syntax:
        --------
        tbldata, outhdr = to_solar(evtdata, hdr)

    """
    

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
    tbldata=evtdata.copy()
    outhdr = hdr.copy()

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



def convert_file(infile, **kwargs):
    """ Wrapper to read the input file, convert the data to solar
        coordinates and then write the output file.
        
        Inputs:
        --------
        Path to the file to be converted. Code assumes this is
        generated by nupipeline and so will have an ".evt" suffix.
        
        Outputs:
        --------
        
        Genreates a new FITS file in the same location as the input
        file but with a _sunpos.evt suffix
        
        Usage:
        --------
        
        convert_file('path/to/file/nuXXXA06_cl.evt')
        

    """
    from astropy.io import fits
    import os, sys
    from os.path import isfile, splitext

    # Make sure input file exists:
    if not isfile(infile):
        print("File does not exist:", infile)
        sys.exit

    # Read in the event data:
    hdulist = fits.open(infile)
    evtdata = hdulist[1].data
    hdr = hdulist[1].header
    hdulist.close()

    # Convert this to solar coordinates
    (newdata, newhdr) = to_solar(evtdata, hdr, **kwargs)

    # # Make the new filename:
    (sfile, ext)=splitext(infile)

    outfile=sfile+'_sunpos.evt'

    # Remove output file if necessary
    if isfile(outfile):
        print(__name__+': '+outfile+' exists! Removing old version...')
        os.remove(outfile)

    print(__name__+': generating file: '+ outfile)
    fits.writeto(outfile, newdata, newhdr)
    return

    
    








