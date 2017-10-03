PRO make_gti, gtifile

IF n_elements(gtifile) EQ 0 THEN gtifile ='gtis.txt'

null = mrdfits('nu20211002001A06_gti.fits', 0, null_header)

gti =  mrdfits('nu20211002001A06_gti.fits', 1, gti_header)

outgti = 'gtis.fits'

mwrfits, null, outgti, null_header, /create

stub = gti[0]
openr, /get_lun, lun, gtifile
WHILE ~(eof(lun)) DO BEGIN
   input = 'string'
   readf, lun, input

   gtis = strsplit(input, ' ' , /extract)

   stub.start=  convert_nustar_time(gtis[0], /from_ut)
   stub.stop = convert_nustar_time(gtis[1], /from_ut)
   
   
   IF n_elements(all_gtis) EQ 0 THEN BEGIN
      all_gtis = stub 
   ENDIF ELSE BEGIN
      all_gtis = [all_gtis, stub]
   ENDELSE


ENDWHILE
close, lun

free_lun, lun
mwrfits, all_gtis, outgti, gti_header





END
