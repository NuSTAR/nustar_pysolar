
fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_B0.dat headfile=header_userbadpixB0.txt outfile=nuBuserbadpix_solar.fits
fmodhead infile=nuBuserbadpix_solar.fits tmpfil=rm_history.txt

fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_B1.dat headfile=header_userbadpixB1.txt outfile=userbadpix_det1.fits
fmodhead infile=userbadpix_det1.fits tmpfil=rm_history.txt
fappend infile=userbadpix_det1.fits+1 outfile=nuBuserbadpix_solar.fits
fmodhead infile=nuBuserbadpix_solar.fits+2 tmpfil=rm_history.txt

fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_B2.dat headfile=header_userbadpixB2.txt outfile=userbadpix_det2.fits
fmodhead infile=userbadpix_det2.fits tmpfil=rm_history.txt
fappend infile=userbadpix_det2.fits+1 outfile=nuBuserbadpix_solar.fits
fmodhead infile=nuBuserbadpix_solar.fits+3 tmpfil=rm_history.txt

fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_B3.dat headfile=header_userbadpixB3.txt outfile=userbadpix_det3.fits
fmodhead infile=userbadpix_det3.fits tmpfil=rm_history.txt
fappend infile=userbadpix_det3.fits+1 outfile=nuBuserbadpix_solar.fits
fmodhead infile=nuBuserbadpix_solar.fits+4 tmpfil=rm_history.txt

fmodhead nuBuserbadpix_solar.fits+0 primary_fpmb.txt
fmodhead nuBuserbadpix_solar.fits+1 header_all.txt
fmodhead nuBuserbadpix_solar.fits+2 header_all.txt
fmodhead nuBuserbadpix_solar.fits+3 header_all.txt
fmodhead nuBuserbadpix_solar.fits+4 header_all.txt
fchecksum nuBuserbadpix_solar.fits update=yes
fverify nuBuserbadpix_solar.fits

rm -f userbadpix_det1.fits userbadpix_det2.fits userbadpix_det3.fits
