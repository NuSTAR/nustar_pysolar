def get_sky_position(time, offset):
    """
    
    Code for converting solar offsets to pointing position.
    
    Synax:
    ----------
    sku_position = get_sky_position(time, offset)

    Input:
    ----------
    time: Date that is parsable by sunpy.time.parse_time()
    
    offset: Offset from the center of the Sun. Must have units from astropy:
    
    i.e.: offset = np.array([1000, 150]) * u.arcsec
    
    
    Output:
    ----------
    sky_position: Two-element array giving the [RA, Dec] coordinates of the 
    
    
    """
    
    # Convert the date into something that's usable by astropy.


    start_date = sunpy.time.parse_time(time)
    astro_time = astropy.time.Time(start_date)
    
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
    
    