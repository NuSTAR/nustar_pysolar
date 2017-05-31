
def convert_nustar_time(t, leap=5):
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
