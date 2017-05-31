from datetime import datetime
from astropy.time import Time

def read_tle_file(tlefile, **kwargs):
    """
    
    Read in a TLE file and return the TLE that is closest to the date you want to
    propagate the orbit to.
    """
    times = []
    line1 = []
    line2 = []
    
    from os import path
    from datetime import datetime
    # Catch if the file can't be opened:
    try:
        f = open(tlefile, 'r')
    except FileNotFoundError:
        print("Unable to open: "+tlefile)
    
    ln=0
    for line in f:
#        print(line)
        if (ln == 0):
            year= int(line[18:20])
            day = int(line[20:23])
            
            times.extend([datetime.strptime("{}:{}".format(year, day), "%y:%j")])
            line1.extend([line.strip()])
            ln=1
        else:
            ln=0
            line2.extend([line.strip()])
    f.close()
    return times, line1, line2

def get_epoch_tle(epoch, tlefile):
    """
    
    Find the TLE that is closest to the epoch you want to search.
    
    epoch is a datetime object, tlefile is the file you want to search through.
    
    """

    times, line1, line2 = read_tle_file(tlefile)
    from datetime import datetime
    from astropy.time import Time
    
    # Allow astropy Time objects
    if type(epoch) is Time:
        epoch = epoch.datetime
        
    mindt = 100.
    min_ind = 0
    for ind, t in enumerate(times):
        dt = abs((epoch -t).days)
        if dt < mindt:
            min_ind = ind
            mindt = dt

    good_line1 = line1[min_ind]
    good_line2 = line2[min_ind]

    return mindt, good_line1, good_line2


   
def convert_nustar_time(t, leap=5):
    ''' 
    
    Converts MET seconds to a datetime object.
    
    Default is to subtract off 5 leap seconds.

    '''
    import astropy.units as u
    mjdref = 55197*u.d
    met = (t - leap)* u.s + mjdref
    met_datetime = Time(met.to(u.d), format = 'mjd').datetime
    return met_datetime


def get_nustar_location(checktime, line1, line2):
    ''' 
    
    Code to determine the spacecraft location from the TLE.
    
    Inputs are a datetime object and the two lines of the TLE you want to use.
    
    Returns a tuple that has the X, Y, and Z geocentric coordinates (in km).
    
    '''
    
    from sgp4.earth_gravity import wgs72
    from sgp4.io import twoline2rv
    from astropy.coordinates import EarthLocation
    
    satellite = twoline2rv(line1, line2, wgs72)
    position, velocity = satellite.propagate(
        checktime.year, checktime.month, checktime.day,
        checktime.hour, checktime.minute, checktime.second)

    return position



def eci2el(x,y,z,dt):
    """
    Convert Earth-Centered Inertial (ECI) cartesian coordinates to ITRS for astropy EarthLocation object.

    Inputs :
    x = ECI X-coordinate 
    y = ECI Y-coordinate 
    z = ECI Z-coordinate 
    dt = UTC time (datetime object)
    """

    from astropy.coordinates import GCRS, ITRS, EarthLocation, CartesianRepresentation
    import astropy.units as u
    
    # convert datetime object to astropy time object
    tt=Time(dt,format='datetime')

    # Read the coordinates in the Geocentric Celestial Reference System
    gcrs = GCRS(CartesianRepresentation(x=x, y=y,z=z), obstime=tt)

    # Convert it to an Earth-fixed frame
    itrs = gcrs.transform_to(ITRS(obstime=tt))

    el = EarthLocation.from_geocentric(itrs.x, itrs.y, itrs.z) 

    return el


    
def get_moon_j2000(epoch, line1, line2, position = None):
    '''
    
    Code to determine the apparent J2000 position for a given
    time and at a given position for the observatory.
    
    epoch needs to be a datetime or Time object.
    
    position is a list/tuple of X/Y/Z positions
    
    '''
    
    from astropy.time import Time
    from astropy.coordinates import get_moon, EarthLocation
    import astropy.units as u
    import sys
    from datetime import datetime
    
    if type(epoch) is Time:
        epoch = epoch.datetime
    
    
    if position is None:
        position = get_nustar_location(epoch, line1, line2)  # position in ECI coords


    t=Time(epoch)
    
    loc = eci2el(*position*u.km,t)

    moon_coords = get_moon(t,loc)
    
    # Get just the coordinates in degrees
    
    ra_moon, dec_moon = moon_coords.ra.degree * u.deg, moon_coords.dec.degree*u.deg
    return ra_moon, dec_moon


