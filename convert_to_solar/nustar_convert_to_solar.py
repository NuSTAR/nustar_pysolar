#!/usr/bin/env python

import sys, getopt, os
from os.path import *

from astropy.io import fits
import astropy.time
import astropy.units as u
from astropy.coordinates import get_sun

from sunpy import sun

import numpy as np

scriptname='nustar_convert_to_solar.py'





def main(argv):
    infile = ''
    tstep = 5.0
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:",["ifile="])
            
    except getopt.GetoptError:
        print(scriptname+'-i <inputfile> -t <timestep in seconds>')
        sys.exit(2)


    for opt, arg in opts:
        if opt == '-h':
            print(scriptname+'-i <inputfile> -t <timestep in seconds>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-t", "--timestep"):
            tstep = float(arg)

# Make the timestep have units of seconds.
    tstep = tstep * u.s
    print(scriptname)
    print("Processing: ", infile)
    print("Updating solar ephemeris every: ",tstep)


### Get the data from the FITS file.
    hdulist = fits.open(infile)

    evtdata = hdulist[1].data
    x = np.array([])
    y = np.array([])

    hdr=hdulist[1].header

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

    x = evtdata['X']
    y = evtdata['Y']
    pi =evtdata['PI']
    met = evtdata['TIME']*u.s

    # Conver the NuSTAR epoch times to astropy datetime objects
    mjdref=hdulist[1].header['MJDREFI']
    time = astropy.time.Time(mjdref*u.d+met, format = 'mjd')

    # Convert X and Y to RA/DEC
    ra_x = ra_ref + (x - x0) * delx / np.cos(dec_ref)
    dec_y = dec_ref + (y - y0) * dely

    print("Loaded: ", len(ra_x), " counts.")
    hdulist.close()


    def get_sun_pos(last_met):
        sun_time = astropy.time.Time(mjdref*u.d+last_met, format = 'mjd')
        astro_sun_pos = get_sun(sun_time)
        
        # Get the center of the Sun, and assign it degrees.
        sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg
        
        # Solar NP roll angle:
        sun_np=sun.solar_north(t=sun_time).cgs
        return sun_pos, sun_np;

# How often you want to update the solar ephemeris:
    tstep = 5. * u.s
    last_met = met[0] - tstep * 2.
    last_i = 0

    sun_x = np.zeros_like(ra_x)
    sun_y = np.zeros_like(dec_y)

    for i in np.arange(len(ra_x)):
        if ( (met[i] - last_met) > tstep ):
            (sun_pos, sun_np) = get_sun_pos(last_met)
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
        offset = sun_pos - ph_pos

    # Account for East->West conversion for +X direction in heliophysics coords
        offset = offset*[-1., 1.]
        
        
        # Project the offset onto the Sun
        delta_offset = ((np.dot(offset, rotMatrix)).to(u.arcsec))
        sun_x[i] = delta_offset[0]
        sun_y[i] = delta_offset[1]

    print("Proccessed: ", i, " of ", len(ra_x))


### Make the new output file

    hdulist = fits.open(infile)

    tbldata=hdulist[1].data
    hdr=hdulist[1].header

    hdulist.close()

    #  change to 0-3000 pixels:
    maxX = 3000
    maxY = 3000
    x0 = maxX / 2.
    y0 = maxY / 2.


    # Header keywords
    for field in hdr.keys():
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
                print(hdr[field][5:8])
                xval = field[5:8]
            if hdr[field] == 'Y':
                print(hdr[field][5:8])
                yval = field[5:8]
    delx = hdr['TCDLT'+xval] * u.deg
    dely = hdr['TCDLT'+yval]*u.deg



    out_sun_x=(sun_x / delx) + x0
    out_sun_y=(sun_y / dely) + y0

    newdelx = delx.to(u.arcsec).value
    newdely = dely.to(u.arcsec).value


    tbldata['X'] = out_sun_x
    tbldata['Y'] = out_sun_y



    hdr['TCRVL'+xval] = 0.
    hdr['TCRPX'+xval] = x0
    hdr['TCDLT'+xval] = delx.to(u.arcsec).value
    hdr['TLMAX'+xval] = maxX
    hdr['TCRVL'+yval] = 0.
    hdr['TCRPX'+yval] = x0
    hdr['TCDLT'+yval] = dely.to(u.arcsec).value
    hdr['TLMAX'+yval] = maxY

    # # Make the new filename:
    (sfile, ext)=splitext(infile)

    outfile=sfile+'_sunpos.evt'

    # Remove output file if necessary
    if isfile(outfile):
        print(outfile, 'exists! Removing old version...')
        os.remove(outfile)

    fits.writeto(outfile, tbldata, hdr)


if __name__ == "__main__":
   main(sys.argv[1:])
