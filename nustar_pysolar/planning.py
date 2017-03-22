from sunpy.time import parse_time
from datetime import timedelta
from datetime import date

from astropy.time import Time
from astropy.coordinates import get_sun
from astropy.coordinates import SkyCoord

from sunpy import sun
import numpy as np
from astropy import units as u



def get_sky_position(time, offset):
    """
    
    Code for converting solar offsets to pointing position.

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
    
    # Convert the date into something that's usable by astropy.


    start_date = parse_time(time)
    astro_time = Time(start_date)
    
    # Use astropy get_sun for Sun sky position.
    # sunpy has a similar function, but it may be giving a different
    # epoch for the RA and dec. We need them in J2000 RA and dec.
    
    astro_sun_pos = get_sun(astro_time)

    # Get the solar north pole angle. cgs --> radians
    sun_np=sun.solar_north(t=time).cgs

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

    return sky_position;
    
def get_nustar_roll(time, angle):
    """
    Code to determine the NuSTAR roll angle for a given field-of-view on the
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
        
    start_date = parse_time(time)
    astro_time = Time(start_date)
    # Get the solar north pole angle. cgs --> radians
    sun_np=sun.solar_north(t=time).deg * u.deg
        
    nustar_roll = np.mod( sun_np + angle, 360*u.deg)

    return nustar_roll;

def download_occultation_times(outdir='./'):
    """
    
    Pull the orbital information from the Berkeley website
    and then returns the location of the file.
     
    Parameters
    ----------
    
    outdir: Optional desired output location. Defaults to the working directory.
    
    Returns
    ----------
    
    Returns the filename that you've downloaded.
    
    Notes
    ---------
    
    Will not work if year_day is in the future. This is the time the file was
    generated, not the date you want to observe. The occultation windows should
    extend roughly a month into the future, though they may be more uncertain
    the further ahead you go.
    
    """
    import os

    # Make sure you've got a trailing slash...
    if not(outdir.endswith('/')):
        outdir+'/'
    
    # Make sure the directory exists and create one if not.
    directory = os.path.dirname(outdir)
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    # Get yesterday's file...
    today = date.today() - timedelta(1) 
    year = str(today.timetuple().tm_year)
    day='{0:03d}'.format(today.timetuple().tm_yday)
    year_doy=year+'_'+day

    
    myname='nustar_pysolar.planning'
    import wget
    
    url='http://hessi.ssl.berkeley.edu/ground_systems/products/NUSTAR/'
    url+=year_doy+'/'
    fname='NUSTAR.'+year_doy+'.SHADOW_ANALYSIS.txt' 
    url+=fname

    # Check to see if the file exists:
    import os.path
    if not(os.path.isfile(outdir+fname)):
        wget.download(url, out=outdir+fname)
    
    return outdir+fname;
    
    
    
    
def _parse_timestamp(tstamp):
    """
    Convenience function for turning the timestamp into a datetime object.
    """
    
    date1 = tstamp.split('/')
    year=date1[0].strip()
    day, time=(date1[1].split())
    
    stub = (year.strip()+'-01-01T00:00:00')
    
    year = parse_time(stub)
    hr, min, sec = time.split(':')
    dt = timedelta(int(day)-1, int(sec), 0, 0, int(min), int(hr))

    return year+dt;


def parse_occultations(infile):
    """
    
    Parse the shadow analysis file to determine the 'in Sun' times. 
    
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

#        print(first, dtfirst)
#        print(second, dtsecond)
#        break
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
    return all_pairs;
    
def sunlight_periods(infile, tstart, tend):
    """
    
    Return the periods when NuSTAR is in Sunlight in the given timerange.
    
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
        
    return in_range;
