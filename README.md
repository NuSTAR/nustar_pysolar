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

### nustar_pysolar

This is a set of python modules that can be used to convert the *NuSTAR* astrophysics data to heliocentric coordinates.

1. Clone the project from github

    $ git clone https://github.com/NuSTAR/nustar_pysolar.git

2. Go back to nustar_pysolar project directory and execute::

    $ python setup.py install


### Documentation

TBD

### notebooks

Contains several jupyter notebooks that demonstrates how to generate a pointing location and *NuSTAR* roll for a given observation. 

The proper order is 

1. (Observation Report)[notebooks/Observation_Report.ipynb]





