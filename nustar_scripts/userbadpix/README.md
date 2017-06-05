# Making the userbadpix files

This contains the scripts to turn off the known bad pixels for *NuSTAR* solar analysis.

You need HEASoft installed to run these scripts (which you should have if you are going
to run the nupipeline).

To check to see if this is install, simply do an:

`echo $HEASDAS`

If you see the path to your heasoft distribution then you're good to go.

If you want to adjust which files are turned off, adjust the lines in each text file.

The columns are:

TSTART     RAWX      RAWY     TYPE

where TSTART should be 78913712 and TYPE should be 4 for all NuSTAR observations.

---

To use these files, give the following to your nupipeline call:

```
fpma_userbpfile=nuAuserbadpix_solar.fits fpmb_userbpfile=nuBuserbadpix_solar.fits
```

...where you've copied the nu[A/B]userbadpix_solar.fits file has been copied to the
directory from which you are going to be calling nupipeline.

