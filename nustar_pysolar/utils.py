
def convert_nustar_time(t, leap=5, astropy_time=False):
    '''Converts MET seconds to a datetime object.
    
    Uses the astropy.time.Time method to do the conversion since you have
    to go to MJD first.
    
    Default is to subtract off 5 leap seconds.

    '''
    import astropy.units as u
    from astropy.time import Time
    
    mjdref = 55197*u.d
    met = (t - leap)* u.s + mjdref
    met_datetime = Time(met.to(u.d), format = 'mjd').datetime
    
    if astropy_time is True:
        from astropy.time import Time
        met_astro = Time(met_datetime)
        return met_astro
    else:
        return met_datetime

def channel_to_keV(pi):
    '''Converts MET seconds to a datetime object.
    
    Uses the astropy.time.Time method to do the conversion since you have
    to go to MJD first.
    
    Default is to subtract off 5 leap seconds.

    '''
    keV = pi * 0.04 + 1.6
    return keV
    
def keV_to_channel(keV):
    '''Converts MET seconds to a datetime object.
    
    Uses the astropy.time.Time method to do the conversion since you have
    to go to MJD first.
    
    Default is to subtract off 5 leap seconds.

    '''
    import numpy as np
    
    pi = np.floor((keV - 1.6) / 0.04)
    
    return pi


def skyfield_ephem(load_path=None, parallax_correction=False, utc=None):
    '''Returns skyfield objects used for doing solar planning and conversion.
    
    Optional keywords:
    
    load_path: Directory location where the bsp files are stored. Defaults to current
    working directory if None.
    
    parallax_corection: Download the latest NuSTAR TLE file and apply the parallax
    correction based on NuSTAR's position in its orbit. Defaults to "False". If "True"
    then uses the nustar_pysolar.io TLE methods to parse the closest TLE entry in the
    NuSTAR database.
    
    If you set parallax_correction=True then 
    
    Returns:
    
    observer, sun, ts
    
    The first two are Skyfield objects. "observer" is either geocentric or the NuSTAR
    location based on the TLE.
    
    "ts" is the Skyfield time object.
    
    '''


    # Initialize Skyfield ephemeris tools.
    from skyfield.api import EarthSatellite, Loader
    from astropy.time import Time
    import sunpy.sun

    if load_path is None:
        load_path = './'
        load=Loader(load_path)
    else:
        load=Loader(load_path)
            
    ts = load.timescale()
    planets = load('de436.bsp')
    earth = planets['Earth']
    sun = planets['Sun']


    if parallax_correction is False:
        observer = earth
    else:
        assert (not utc is None),"Must set UTC when using parallax correction!"            
        import nustar_pysolar.io as io
        tlefile = io.download_tle(outdir=load_path)
        mindt, line1, line2 = io.get_epoch_tle(utc, tlefile)
        nustar = EarthSatellite(line1, line2)
        observer = earth + nustar

    ts = load.timescale()


    return observer, sun, ts