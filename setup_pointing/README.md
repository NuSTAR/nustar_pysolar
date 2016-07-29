
## Overview

This scripts gives an example of how to generate pointing positons for *NuSTAR* (in RA and dec) based on pointing position offsets on the solar disk.

### Libraries we'll use

Here we require numpy, astropy, and sunpy.

sunpy is mostly used for date convenience and to get the Sun's north pole angle with respect to celestial north.


```python
from sunpy import sun
import sunpy.time

import astropy.time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import get_sun

import numpy as np
```

### Sky Position Code
The code below defines the function to convert a given offset on the solar disk into an RA/Dec pointing position.

The code assumes that time is a datetime object or a string that looks like '2016-07-03T10:33:10' and that offset is a two element numpy vector like [1000, 150] where the values are the X and Y offsets from the center of the Sun.

Here we also assume that offset has a unit assocaited with it (i.e. u.arcsec) so that we don't have to worry about hardcording any unit conversions. Just in case we ever want to do arcminute offsets or the like.


```python
def get_sky_position(time, offset):
    "Code for converting solar offsets to pointing position."

    

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
```

### NuSTAR Roll angle code

This code just computes the "Sky" PA angle required for NuSTAR to have a given roll.

Some common roll angles you'll want to use:


**Square** DET0 at NE / SE / SW / NW of FoV, angle = 0 / 90 / 180 / 270 degrees.

**Diamond** (Det 0 at the soouth): angle = 45






```python
def get_nustar_roll(time, angle):
    "Code to get NuSTAR roll for given FoV."
        
    start_date = sunpy.time.parse_time(time)
    astro_time = astropy.time.Time(start_date)
    # Get the solar north pole angle. cgs --> radians
    sun_np=sun.solar_north(t=time).deg * u.deg
        
    nustar_roll = np.mod( sun_np + angle, 360*u.deg)

    return nustar_roll;
```

### Usage
The time input should be a date string or a datetime object.


```python
aim_time='2016-07-26T19:53:15.00'
offset = np.array([1000, 150]) * u.arcsec
sky_pos = get_sky_position(aim_time, offset)
angle=90 * u.deg
nustar_roll = get_nustar_roll(aim_time, angle)



print("Aim time: ", aim_time)
print("Coordinates: ",sky_pos)
print("NuSTAR Sky PA: ",nustar_roll)



```

    Aim time:  2016-07-26T19:53:15.00
    Coordinates:  [ 126.04053869   19.33666449] deg
    NuSTAR Sky PA:  98.8658844320891 deg


### Slightly more advanced usage:

Here if you want to specify the dwell start/stop times and compute the aim time automatically.


```python
dwell_start='2016-07-26T19:22:10'
dwell_end='2016-07-26T20:24:20'

start_date = sunpy.time.parse_time(dwell_start)
end_date = sunpy.time.parse_time(dwell_end)

dt = end_date-start_date
aim_time = start_date + dt*0.5


offset = np.array([1000, 150]) * u.arcsec
sky_pos = get_sky_position(aim_time, offset)
angle=90 * u.deg
nustar_roll = get_nustar_roll(aim_time, angle)


print("Aim time: ", aim_time)
print("Coordinates", sky_pos)
print("NuSTAR Sky PA: ",nustar_roll)

```

    Aim time:  2016-07-26 19:53:15
    Coordinates [ 126.04053869   19.33666449] deg
    NuSTAR Sky PA:  98.8658844320891 deg

