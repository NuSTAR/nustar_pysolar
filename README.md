# *NuSTAR* Solar Python Repo

## Overview:

**This is currently work in progress. If you encounter issues, let us know.**

This repo contains scripts for converting *NuSTAR* solar observations to usable
heliocentric coordinates.

This branch uses python-based code rather than IDL scripts so that we can take advantage
of the wonderful open source developments in astropy and sunpy.

All python scripts shown here are written for python 3.5.

Library Requirements (at minimum, and with their dependencies):
  astopy
  numpy
  sunpy

We recommend using [Anaconda](https://www.continuum.io/downloads) for installation
of astropy/numpy/everything else via conda.

See the [sunpy documenation](http://sunpy.org) for details of how to install sunpy
via conda.

See Iain Hannah's overview figures for the solar observations (made using SSWIDL)
[here](http://ianan.github.io/nsigh_all/).

Interested in helping out or adding code? Feedback is great, report any problems via 
`issues`_ page.

If you have code you might want to contribute, get in touch with me to join the Slack
group and/or issue a pull request.

## Contents: 

### nustar_pysolar

This is a set of python modules that can be used to convert the *NuSTAR* astrophysics
data to heliocentric coordinates.

1. Clone the project from github:

>    git clone https://github.com/NuSTAR/nustar_pysolar.git

2. Go to the nustar_pysolar project directory and execute:

>    python setup.py install


### Documentation

To build the documentation do the following:

1. cd to the directory where you have cloned nustar_pysolar

2. Issue the command:

>   build-sphinx docs/ docs/_build

The documentation will then be built into the nustar_pysolar/docs/\_build directory. The
top page is index.html.


### notebooks

We've provided several jupyter notebooks that demonstrates how to convert the *NuSTAR*
solar data into heliocentric coordinates and data formats that align with Sunpu.

This series of notebooks should result in you have a sunpy-style map object that can
then be combined with other sunpy objects (AIA, etc) for making pretty pictures.

The proper order is:

1. [Observation Report](notebooks/Observation_Report.ipynb)
2. [Convert Example](notebooks/Convert_Example.ipynb)
3. [Map Example](notebooks/Map_Example.ipynb)

## How to get *NuSTAR* solar data

You can search the [*NuSTAR* table](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dnumaster&Action=More+Options) at the HEASARC, using SOL in the obs_type search bar to locate all of the *NuSTAR* solar observations. Any data that are public can be downloaded and reprocessed for solar work.

If you have data that you'd like to see that aren't public let, let me know.





