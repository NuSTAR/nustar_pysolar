from datetime import timedelta

from astropy import units as u
import numpy as np
from sunpy.time import parse_time


def get_sky_position(time, offset):
    """Code for converting solar offsets to pointing position.

    Parameters
    ----------
    time: Date that is parsable by sunpy.time.parse_time()

    i.e.,

    time='2016-07-26T19:53:15.00'

    offset: Offset from the center of the Sun. Must have units from astropy:

    i.e.: offset = np.array([1000, 150]) * u.arcsec


    Returns
    ----------
    sky_position: Two-element array giving the [RA, Dec] coordinates of the


    Notes
    ----------
    Syntax:

    sky_position = get_sky_position(time, offset)



    """

    from astropy.coordinates import get_sun
    from astropy.time import Time
#     Replaced with newer sunpy v1 function
#     from sunpy import sun
    from sunpy.coordinates import sun

    # Convert the date into something that's usable by astropy.


    start_date = parse_time(time)
    astro_time = Time(start_date)

    # Use astropy get_sun for Sun sky position.
    # sunpy has a similar function, but it may be giving a different
    # epoch for the RA and dec. We need them in J2000 RA and dec.

    astro_sun_pos = get_sun(astro_time)

    # Get the solar north pole angle. cgs --> radians
#     Update for sunpy v1.0+
#     sun_np=sun.solar_north(t=time).cgs
    sun_np=sun.P(time).cgs

    # Get the center of the Sun, and assign it degrees.
    # Doing it this was is necessary to do the vector math below.
    sun_pos = np.array([astro_sun_pos.ra.deg, astro_sun_pos.dec.deg])* u.deg

    # Rotation matrix for a counter-clockwise rotation since we're going
    # back to celestial north from solar north
    rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
                         [-np.sin(sun_np),  np.cos(sun_np)]])

    # Project the offset onto the Sun
    delta_offset = np.dot(offset, rotMatrix)

    # Scale to RA based on the declination.
    delta_offset = delta_offset * np.array([1. / np.cos(sun_pos[1]), 1.])

    # Account for the fact that +Ra == East and we have defined +X = West
    delta_offset = delta_offset * [-1.0, 1.0]

    # Apply the offset and return the sky position.
    sky_position = sun_pos + delta_offset

    return sky_position


def get_skyfield_position(time, offset, load_path=None, parallax_correction=False):
    """Code for converting solar coordinates to astrometric (J200) RA/Dec coordinates.

    Parameters
    ----------
    time: Date that is parsable by sunpy.time.parse_time()

    i.e.,

    time='2016-07-26T19:53:15.00'

    offset: Offset from the center of the Sun. Must have units from astropy:

    i.e.: offset = np.array([1000, 150]) * u.arcsec


    load_path (optional): Relative path from currently location to store bsp files

    parallax_correction: Use the NuSTAR TLE to correct for orbital parallax


    Returns
    ----------
    sky_position: Two-element array giving the [RA, Dec] coordinates of the
    target location. Note this is given in astrometric (J2000) RA/Dec, which is what
    we need for the NuSTAR planning system.

    Notes
    ----------
    Syntax:

    skyfield_position = get_skyfield_position(time, offset)

    """

    from astropy.time import Time
#     Replaced with newer sunpy v1 function
#     from sunpy import sun
    from sunpy.coordinates import sun
    from nustar_pysolar.utils import skyfield_ephem
    start_date = parse_time(time)
    utc = Time(start_date)

    observer, sunephem, ts = skyfield_ephem(load_path=load_path,
                                        parallax_correction=parallax_correction,
                                        utc=utc)

    tcheck = ts.from_astropy(utc)
    geocentric = observer.at(tcheck).observe(sunephem)
    this_ra_geo, this_dec_geo, dist = geocentric.radec()


    # Get the solar north pole angle. cgs --> radians
#     sun_np = sunpy.sun.solar_north(t=time).cgs
    #     Update for sunpy v1.0+
    sun_np=sun.P(time).cgs

    # Get the center of the Sun, and assign it degrees.
    # Doing it this was is necessary to do the vector math below.
    sun_pos = np.array([this_ra_geo.to(u.deg).value, this_dec_geo.to(u.deg).value])*u.deg


    # Rotation matrix for a counter-clockwise rotation since we're going
    # back to celestial north from solar north
    rotMatrix = np.array([[np.cos(sun_np), np.sin(sun_np)],
                          [-np.sin(sun_np), np.cos(sun_np)]])

    # Project the offset onto the Sun
    delta_offset = np.dot(offset, rotMatrix)

    # Scale to RA based on the declination.
    delta_offset = delta_offset * np.array([1. / np.cos(sun_pos[1]), 1.])

    # Account for the fact that +Ra == East and we have defined +X = West
    delta_offset = delta_offset * [-1.0, 1.0]

    # Apply the offset and return the sky position.
    sky_position = sun_pos + delta_offset

    return sky_position


def get_nustar_roll(time, angle):
    """Code to determine the NuSTAR roll angle for a given field-of-view on the
    Sun for a given time.
    
    Parameters
    ----------
    time: Date that is parsable by sunpy.time.parse_time()
    
    i.e. 
    
    time='2016-07-26T19:53:15.00'
    
    
    angle: Desired roll offset from solar north in degrees.
    
    For a "square" field of view, use angle=0 / 90 / 180 / 270 to have DET0
    at the NE / SE / SW / NW corners of a square field of view.
    
    For a "diamond" with DET0 to the south, use angle = 45.
        
    Returns
    ----------
    nustar_roll: NuSTAR PA angle with respect to celestial north.
    
    """

#     Replaced with newer sunpy v1 function
#     from sunpy import sun
    from sunpy.coordinates import sun

    # Get the solar north pole angle. cgs --> radians
#     sun_np=sun.solar_north(t=time).deg * u.deg
    #     Update for sunpy v1.0+
    sun_np=sun.P(time).deg*u.deg

    nustar_roll = np.mod(sun_np + angle, 360*u.deg)

    return nustar_roll;

def _parse_timestamp(tstamp):
    """Convenience function for turning the SOC timestamp into a datetime object.
    """
    
    date1 = tstamp.split('/')
    year=date1[0].strip()
    day, time=(date1[1].split())
    
    stub = (year.strip()+'-01-01T00:00:00')
    
    year = parse_time(stub)
    hr, min, sec = time.split(':')
    dt = timedelta(int(day)-1, int(sec), 0, 0, int(min), int(hr))

    return year+dt;
    
def _parse_SOC_timestamp(tstamp):
    """Convenience function for turning the timestamp into a datetime object.
    """
    
    date1 = tstamp.split(':')

    year = date1[0]
    day = date1[1]
    hr = date1[2]
    min = date1[3]
    sec = date1[4]
    
    
    
    stub = (year.strip()+'-01-01T00:00:00')
    
    year = parse_time(stub)
#    hr, min, sec = date1[2:4]
    dt = timedelta(int(day)-1, int(sec), 0, 0, int(min), int(hr))
    return year+dt;

def parse_occultations(infile):
    """Parse the shadow analysis file to determine the 'in Sun' times. 
    
    Parameters
    ----------
    
    infile: Input file to be parsed.
    
    
    Returns
    ----------
    
    Returns a list of [ [start, stop], [start stop] ] times where start means
    you egress from Earth shadow into the sunlight, while stop means you
    re-enter Earth shadow.
    
    Notes
    ---------

    
    """

    f = open(infile)
    all_pairs = []
    start = 0
    for ind,line in enumerate(f):

        # Little parser here to find the right place to start reading in...
        if (line.find("Shadow Begin") != -1):
            start=start+1

        # Skips over additional lines of whitespace.
        if(start == 0):
            continue
        if(start <3):
            start+=1
            continue
        
        # Get the first date string:
        fields = line.split('-')
        
        first = fields[0]
        dtfirst = _parse_timestamp(first)
    
    
        second = (fields[1].split('UTC'))[0].strip()
        dtsecond=_parse_timestamp(second)

        # Since the file actually gives the start/stop times of going into
        # earthshadow, we actually want the "In Sun" times, which is the egress
        # from earthshadow and the entry into the next earthshadow.

        # Note that this skips the first row.
        if(start == 3):
            start+=1
        else:
            all_pairs.append([last, dtfirst])
            
        # Store the last entry to add in the next time around...
        last=dtsecond

    f.close()
    return all_pairs

def sunlight_periods(infile, tstart, tend):
    """Return the periods when NuSTAR is in Sunlight in the given timerange.
    
    Parameters
    ----------
    
    tstart, tend: ISO formatted times or something else that
    sunpy.time.parse_time() can read.
    
    i.e.
    
    tstart='2017-03-11T23:09:10'
        
    infile: Input file to be parsed. This should the value returned by
    nustar_pysolar.download_occultation_times()
        
    Returns
    ----------
    
    Returns a list of [ [start, stop], [start stop] ] times where start means
    you egress from Earth shadow into the sunlight, while stop means you
    re-enter Earth shadow.
    
    The list has been filtered to only include those epochs that span the given
    time range.
    
    Notes
    ---------

    """
    import os.path
    if not(os.path.isfile(infile)):
        print('Error in nustar_pysolar.sunlight_periods.')
        print('Input file: '+infile+' does not exist.')
        return -1;

    
    all_pairs = parse_occultations(infile)
    checkstart = parse_time(tstart)
    checkend = parse_time(tend)
    in_range = []

    
    set=0
    for pair in all_pairs:
        dtmin = (pair[0] - checkstart)
        dtmax = (pair[1] - checkstart)
        if ( (pair[1] > checkstart) ):
            set=1
        if (set == 0):
            continue
        if ( pair[1] > checkend ):
            break
        in_range.append(pair)
        
    if len(in_range) == 0:
        print('Error in function: '+sunlight_periods.__name__)
        print('No dates found in range. Pick a different occultation file.')
        return -1
    else:
        return in_range


def make_mosaic(orbit, outfile='mosaic.txt', write_output=False, make_regions=False,
                reg_pref='testbox', extra_roll=0.*u.deg, write_sun=False):
    ''' 
    Code to make a mosaic for a 5x5 tiled array on the Sun.
    
    Input:
    
    tstart = '2018-05-28T15:37:00'
    tend = '2018-05-28T23:10:00'
    
    positions = make_mosaic(tstart, tend, write_output=True)
    
    Optional flags:
    
    write_output = [False] / True
        Write the output pointing positions in NuSTAR SOC readable formats in 'outfile' for all of the pointings.

    outfile = ['mosaic.txt']
        Output file if write_output is used.
    
    make_regions: [False] / True
        Make ds9 region files for each tile so that you can see how the FoV moves with each mosaic location.

    reg_pref: 'testbox'
        The prefix for the region files. Useful if you want to make this meaningful.

    Output mosaic file has columns of:
    "Arrive By Time"     RA      DEC   RA_SUN  DEC_SUN

    '''
    import numpy as np
    box_pa = get_nustar_roll(orbit[0], extra_roll)
    sun_pa = get_nustar_roll(orbit[0], 0.)
    pa = box_pa + 90*u.deg
    

    base = np.array([-1.45, -0.725, 0, 0.725, 1.45])
    xsteps = np.append(base, np.flip(base, 0))
    xsteps = np.append(xsteps, base)
    xsteps = np.append(xsteps, np.flip(base, 0))
    xsteps = np.append(xsteps, base)

    ysteps = np.array(np.zeros(5) + 1.45)
    ysteps = np.append(ysteps, np.zeros(5) + 0.725)
    ysteps = np.append(ysteps, np.zeros(5))
    ysteps = np.append(ysteps, np.zeros(5)-0.725)
    ysteps = np.append(ysteps, np.zeros(5)-1.45)


    # Rotation matrix for a clockwise rotation on the solar disk:
    rotMatrix = np.array([[np.cos(extra_roll), np.sin(extra_roll)],
                         [-np.sin(extra_roll),  np.cos(extra_roll)]])

    
    dt = (orbit[1] - orbit[0]) / 25.

    print("Orbit start: {} Orbit end: {}".format(orbit[0].isoformat(), orbit[1].isoformat()))
    print("Dwell per position:", dt.total_seconds())
    print("")
    print("NuSTAR Roll Angle to get roll relative to Sun of {:.02f} is {:.02f} deg".format(extra_roll.value, box_pa.value))
    print("Step of FOV PA direction is {:.02f} deg".format(pa.value))
    print("")


    if write_output is True:
        f = open(outfile, 'w')
    
    aim_time = orbit[0]
    for ind, pair in enumerate(zip(xsteps, ysteps)):
        arrive_time = aim_time
        aim_time = aim_time + dt
        
        
        # Make this 10-arcmin steps.
        step_size = 10 * u.arcmin

        
        # Sun-center location:
        offset = [0., 0.]*u.deg
        sun_pos = get_skyfield_position(aim_time, offset, load_path='../data', parallax_correction=True)
#        print('Sun time: {} RA (deg): {} Dec (deg): {}'.format(aim_time.isoformat(), sun_pos[0], sun_pos[1]))
        
        # Pointing location        
        # Rotate to the correct orientation on the solar disk:
        offset = (np.dot(pair, rotMatrix) * step_size).to(u.deg)        
        sky_pos = get_skyfield_position(aim_time, offset, load_path='../data', parallax_correction=True)
#        print('Arrive time: {} RA (deg): {} Dec (deg): {}'.format(arrive_time.isoformat(), sky_pos[0], sky_pos[1]))

        if make_regions:
            make_test_region(sky_pos[0], sky_pos[1], box_pa, sun_pos[0], sun_pos[1],sun_pa, outname=reg_pref+'{}.reg'.format(ind))
        if write_output:
            if write_sun:
                f.write('{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f}\n'.format(arrive_time.strftime('%Y:%j:%H:%M:%S'),
                                                       sky_pos[0].value, sky_pos[1].value,
                                                        sun_pos[0].value, sun_pos[1].value))
            else:
                f.write('{0} {1:.4f} {2:.4f} \n'.format(arrive_time.strftime('%Y:%j:%H:%M:%S'),
                                       sky_pos[0].value, sky_pos[1].value))


    if write_output:
        f.close()


def make_test_region(boxra, boxdec, boxpa, sunra, sundec, sunpa, outname='testbox.reg'):
    ''' 
    Code to produce a ds9 region file showing the Sun and the FoV on the sky.
    
    Inputs:
    
    boxra, boxdec: These are the center of the FoV. Must have astropy units.
    boxpa: This is the SKY PA angle. Must have astropy units.
    
    sunra, sundec: The center of the Sun.
    sunpa: The solar north pole angle
    
    outname: the name of the region file that you want to produce.
    Defaults to 'tetsbox.reg'

    The output is a ds9 region file. Go into ds9, enter the Sun ra/dec coordinates and
    load an optical stellar background field (under the 'Analysis' tab), then load
    the region file to see what's going on.

    '''   
    f = open(outname, 'w')
    f.write("# Region file format: DS9 version 4.1\n")
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 \
source=1 \n')
    f.write('fk5 \n')
    
    outstring = 'box({0}, {1}, 720", 720", {2})'.format(boxra.value, boxdec.value, boxpa.value+90 % 360)
    f.write(outstring+'\n')
    
    outstring = 'circle({0}, {1}, 960.5")'.format(sunra.value, sundec.value)


    f.write(outstring+'\n')

    outstring = 'vector({0}, {1}, 960.5", {2})'.format(sunra.value, sundec.value, sunpa.value+90 % 360 )
    f.write(outstring+'\n')

    f.close
    return


