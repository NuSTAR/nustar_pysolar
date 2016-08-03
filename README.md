# *NuSTAR* Solar Repo

## Overview:

**This is currently work in progress.**

This repo contains scripts for converting *NuSTAR* solar observations to usable heliocentric coordinates.

This branch uses python-based code rather than IDL scripts so that we can take advantage of the wonderful open source developments in astropy and sunpy.

All python scripts shown here are written for python 3.5.

Library Requirements (at minimum, and with their dependencies):
  astopy
  numpy
  sunpy

We recommend using [Anaconda](https://www.continuum.io/downloads) for installation of astropy/numpy/everything else via conda.

See the [sunpy documenation](http://sunpy.org) for details of how to install sunpy via conda.

See Iain Hannah's overview figures for the solar observations (made using SSWIDL) [here](http://ianan.github.io/nsigh_all/).

## Contents: 

### setup_pointing:

Contains a jupyter notebook that demonstrates how to generate a pointing location and *NuSTAR* roll for a given observation. 

### convert_to_solar:

Contains code to convert the output of nupipeline into heliophysics coordinates. See the jupyter notebook for documentation and an example. Also has a standalone python script to convert individual files.

This directory also contains a jupyter notebook and associated standalone script to produce sunpy-compatible map files.

### solar_analysis:

Empty. Will probably contain sunpy specific (or [SolarSoft](http://www.lmsal.com/solarsoft/) equivalent) example code to make pretty figures, etc.






