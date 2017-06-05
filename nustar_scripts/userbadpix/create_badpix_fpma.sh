
fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_A0.dat headfile=header_userbadpixA0.txt outfile=nuAuserbadpix_solar.fits
fmodhead infile=nuAuserbadpix_solar.fits tmpfil=rm_history.txt

fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_A1.dat headfile=header_userbadpixA1.txt outfile=userbadpix_det1.fits
fmodhead infile=userbadpix_det1.fits tmpfil=rm_history.txt
fappend infile=userbadpix_det1.fits+1 outfile=nuAuserbadpix_solar.fits
fmodhead infile=nuAuserbadpix_solar.fits+2 tmpfil=rm_history.txt

fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_A2.dat headfile=header_userbadpixA2.txt outfile=userbadpix_det2.fits
fmodhead infile=userbadpix_det2.fits tmpfil=rm_history.txt
fappend infile=userbadpix_det2.fits+1 outfile=nuAuserbadpix_solar.fits
fmodhead infile=nuAuserbadpix_solar.fits+3 tmpfil=rm_history.txt

fcreate cdfile=columns_userbadpix.txt datafile=data_userbadpix_A3.dat headfile=header_userbadpixA3.txt outfile=userbadpix_det3.fits
fmodhead infile=userbadpix_det3.fits tmpfil=rm_history.txt
fappend infile=userbadpix_det3.fits+1 outfile=nuAuserbadpix_solar.fits
fmodhead infile=nuAuserbadpix_solar.fits+4 tmpfil=rm_history.txt

fmodhead nuAuserbadpix_solar.fits+0 primary_fpma.txt
fmodhead nuAuserbadpix_solar.fits+1 header_all.txt
fmodhead nuAuserbadpix_solar.fits+2 header_all.txt
fmodhead nuAuserbadpix_solar.fits+3 header_all.txt
fmodhead nuAuserbadpix_solar.fits+4 header_all.txt
fchecksum nuAuserbadpix_solar.fits update=yes
fverify nuAuserbadpix_solar.fits

rm -f userbadpix_det1.fits userbadpix_det2.fits userbadpix_det3.fits
