This is a temporary IDL placeholder until we re-write this into python. You're required to have some version of the
AstroLib since this uses the "mrdfits" and "mwrfits" IDL wrappers. See the instructions in the
[nustar-idl](https://github.com/NuSTAR/nustar-idl) repository if you don't already have that installed.

"gti.txt" is the input time (in standard UTC timecodes).

"gti.fits" is the output

"nu20211002001A06_gti.fits" is a template GTI file that I'm using to have the correct header keywords.

"make_gtis.pro" is the IDL script that you should run to produce the "gti.fits" output. It relies on the
time conversion script "convert_nustar_time.pro" to change the UT time codes into NuSTAR MET times. Again,
this is included as part of the nustar-idl repo, but I've included it here for expendiency.

You then feed the "gti.fits" file to all NuSTAR scripts via the "usrgtifile" keyword. See the
[nuproducts helpfile](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/caldb/help/nuproducts.html)
and/or the [nupipeline helpfile](https://heasarc.gsfc.nasa.gov/ftools/caldb/help/nupipeline.html) for further
instructions on how to use the usrgti files.
