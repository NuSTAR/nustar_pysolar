NuSTAR Solar Python Utilities
=============================

|Readthedocs| |Astropy|


Overview:
--------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nuclear Spectroscopic Telescope ARray
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://www.nustar.caltech.edu/system/avm_image_sqls/binaries/26/page/nustar_artistconcept_2.jpg?1393022433
    :target: http://www.nustar.caltech.edu
    :alt: NuSTAR

Installation
------------

This repo requires some flavor of Python 3.

Base installation is accomplished via pip, but we recommend installing into its own
conda environment first.


.. code-block:: bash

    conda config --add channels conda-forge 
    conda config --set channel_priority strict
    conda create --name nustar_pysolar python=3.9
    git clone https://github.com/NuSTAR/nustar_pysolar.git


Then install the development version of this using pip. This will also install all
required dependencies. Move to where you downloaded nustar_pysolar, and then do this:

.. code-block:: bash

    pip install -e .

Solar Observations
-------------------

See Iain Hannah's overview figures for the solar observations (made using SSWIDL)
`here <http://ianan.github.io/nsigh_all/>`_.

Interested in helping out or adding code? Feedback is great.

If you have code you might want to contribute, get in touch with me to join the Slack
group and/or issue a pull request.


Notebooks
----------

We've provided several jupyter notebooks that demonstrates how to convert the *NuSTAR*
solar data into heliocentric coordinates and data formats that align with Sunpy.


This series of notebooks should result in you have a sunpy-style map object that can
then be combined with other sunpy objects (AIA, etc) for making pretty pictures.

The proper order is:

1. `Observation Report <notebooks/Observation_Report.ipynb>`_
2. `Convert Example <notebooks/Convert_Example.ipynb>`_
3. `Map Example <notebooks/Map_Example.ipynb>`_

We also have provided an example notebook for how to plan a *NuSTAR* solar observation

1. `Planning example <notebooks/Planning_Example.ipynb>`_

How to get *NuSTAR* solar data
-------------------------------


You can search the `NuMASTER table <https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dnumaster&Action=More+Options>`_ at the HEASARC, using SOL in the obs_type search bar to locate all of the *NuSTAR* solar observations. Any data that are public can be downloaded and reprocessed for solar work.




.. image:: https://readthedocs.org/projects/nustar-pysolar/badge/?version=latest
    :target: https://nustar-pysolar.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
.. |Astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge
