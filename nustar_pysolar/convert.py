import astropy.time
import astropy.units as u
from astropy.coordinates import get_sun
# don't think the following is needed
# import sunpy.map
# don't think the following is needed
# from sunpy import sun
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
    # Deprecated from sunpy v1 onwards
#     from sunpy.coordinates import get_sun_P
#  Use new version instead
    from sunpy.coordinates import sun

    sun_time = astropy.time.Time(met, format = 'mjd')
    astro_sun_pos = get_sun(sun_time)

    # Get the center of the Sun, and assign it degrees.
    sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg

    # Solar NP roll angle:
#     sun_np = get_sun_P(last_met)
    sun_np = sun.P(met)

    return sun_pos, sun_np;



def _xy_to_radec(evtdata, hdr):
    """ Conversion function to go from X/Y coordinates
        in the FITS file to RA/Dec coordinates.
    """
    
    from nustar_pysolar.utils import convert_nustar_time


# Parse the header information
    for field in hdr.keys():
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
#                print(hdr[field][5:8])
                xval = field[5:8]
            if hdr[field] == 'Y':
#                print(hdr[field][5:8])
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
    met = convert_nustar_time(evtdata['Time'], astropy_time=True)
    
#     mjdref=hdr['MJDREFI']
#     met = evtdata['TIME']*u.s + mjdref*u.d

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


def _delta_solar_skyfield(ra_x, dec_y, met, **kwargs):
    """ Function to compute the offsets from the center of
        the Sun as a function of time.

        Use the tStep argument to define how often you want
        to update the solar ephemeris. Default is every 5 seconds.

        Inputs: ra_x, dec_y, and met are all arrays that contain the
        RA, Dec, and time of arrival of each count, respectively. The arrival time
        must have astropy units attached to it.

        Outputs: sun_x, sun_y are the x and y values (in arcseconds)
        from the center of the Sun in the +North and +West directions.

    """
    import astropy.units as u
    from nustar_pysolar.utils import skyfield_ephem
#      Don't think this is needed
#     from sunpy import sun
    # Deprecated from sunpy v1 onwards
#     from sunpy.coordinates import get_sun_P
#  Use new version instead
    from sunpy.coordinates import sun
    
    
    # How often you want to update the solar ephemeris:
    tStep=kwargs.get('tStep', 5.0)
    tStep = tStep * u.s

    load_path=kwargs.get('load_path', None)
    observer, TheSun, ts = skyfield_ephem(load_path = load_path, 
                                        parallax_correction=True,
                                        utc=met[0])
    
        

    # How many events do you want to do?
    maxEvt=kwargs.get('maxEvt', len(ra_x))

    # Keep last time we updated things
    last_met = met[0] - tStep * 2.
    last_i = 0
    
    sun_x = np.zeros_like(ra_x)
    sun_y = np.zeros_like(dec_y)
    for i in np.arange(len(ra_x)):
        if( (met[i] - last_met) > tStep ):
            last_met = met[i]
            

            tcheck = ts.from_astropy(last_met)
            astrometric = observer.at(tcheck).observe(TheSun)
            this_ra, this_dec, dist = astrometric.radec()
            
            

            # Get the center of the Sun, and assign it degrees.
            # Doing it this was is necessary to do the vector math below.
            sun_pos = np.array([this_ra.to(u.deg).value,
                this_dec.to(u.deg).value])*u.deg

#             sun_np = get_sun_P(last_met)
            sun_np = sun.P(last_met)

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
    import astropy.time
    import astropy.units as u
    from astropy.coordinates import get_sun
    # don't think the following is needed
    # import sunpy.map
    # don't think the following is needed
    # from sunpy import sun
    import numpy as np
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
    (sun_x, sun_y) = _delta_solar_skyfield(ra_x, dec_y, met, **kwargs)


    # Parse the header information to get the native bin size
    for field in list(hdr.keys()):
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
#                print(hdr[field][5:8])
                xval = field[5:8]
            if hdr[field] == 'Y':
#                print(hdr[field][5:8])
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
    #     tbldata['X'] = out_sun_x                               **Removed
    #     tbldata['Y'] = out_sun_y                               **Removed

    #     #Update header information                             **Removed
    #     outhdr['TCRVL'+xval] = 0.                              **Removed
    #     outhdr['TCRPX'+xval] = x0                              **Removed
    #     outhdr['TCDLT'+xval] = 1.0 * delx.to(u.arcsec).value   **Removed
    #     outhdr['TLMAX'+xval] = maxX                            **Removed
    #     outhdr['TCRVL'+yval] = 0.                              **Removed
    #     outhdr['TCRPX'+yval] = x0                              **Removed
    #     outhdr['TCDLT'+yval] = dely.to(u.arcsec).value         **Removed
    #     outhdr['TLMAX'+yval] = maxY                            **Removed
    
    #     return tbldata, outhdr                                 **Removed
    
    #**************************** added ***************************************************
    
    # Astropy update (>3. I think) means that the keywords
    # that we need to change are now protected.
    # This means that they are now determined from the data
    # columns and can't just be reassigned.
    # **See: https://github.com/astropy/astropy/issues/7145
    # This is why lines 70--78 were removed.
    
    from astropy.io import fits
    # Remove the keywords so they don't have to be saved out to be updated.
    del outhdr['TCRVL'+xval]
    del outhdr['TCRPX'+xval]
    del outhdr['TCDLT'+xval] 
    
    del outhdr['TCRVL'+yval] 
    del outhdr['TCRPX'+yval] 
    del outhdr['TCDLT'+yval] 
    
    # These ones can still be changed the old way for some reason.
    outhdr['TLMAX'+xval] = maxX
    outhdr['TLMAX'+yval] = maxY
    
    # Get columns from the data
    orig_cols = tbldata.columns
    
    # Make new columns for the solar position coordinates in the X and Y fields.
    # This method won't overwrite the X and Y fields if they are already there. 
    # This is why line 67 and 68 were removed.
    # To create the columns, the fooloowing format is used.
    # **See: https://docs.astropy.org/en/stable/io/fits/usage/table.html
    # fits.Column(name=TTYPE, format=TFORM, unit=TUNIT, null=TNULL, 
    #             coord_ref_point=TCRPX, coord_ref_value=TCRVL, coord_inc=TCDLT, 
    #             coord_type=TCTYP, coord_unit=TCUNI,
    #             array=data)
    
    # Remove the RA---DEC X and Y data fields from the original columns
    orig_cols.del_col('X')
    orig_cols.del_col('Y')
    
    # Now create the new columns which contain the keywords we
    # need as well as the solar X and Y coordinates
    new_cols = fits.ColDefs([fits.Column(name='X', format='1I', unit='pixel', null=-1, 
                                         coord_ref_point=x0, coord_ref_value=0., 
                                         coord_inc=1.0 * delx.to(u.arcsec).value, 
                                         coord_type="Heliopro", coord_unit="arcsec",
                                         array=out_sun_x),
                             fits.Column(name='Y', format='1I', unit='pixel', null=-1, 
                                         coord_ref_point=y0, coord_ref_value=0., 
                                         coord_inc=dely.to(u.arcsec).value, 
                                         coord_type="Heliopro", coord_unit="arcsec",
                                         array=out_sun_y)])
    
    # Combine the old and new columns to create a new hdu structure
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    
    # separate the new hdu intp its data and header components
    # to be able to be consistent with to_solar()
    new_tbldata, hdr_from_new_columns = hdu.data, hdu.header
    
    # return the new table data, but also combine the original header with the new updated one
    return new_tbldata, outhdr + hdr_from_new_columns









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

    
    








