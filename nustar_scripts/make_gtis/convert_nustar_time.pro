FUNCTION convert_nustar_time, time, $
                              from_ut = from_ut, $
                              from_mjd = from_mjd , $
                              from_nustar =from_nustar, $
                              from_aft=from_aft, $
                              nustar=nustar, $
                              mjd = mjd, $
                              jd=jd,$
                              fits=fits, $
                              ut = ut, $
                              aft = aft


; Routine for converting times to/from NuSTAR epoch times.
; Inputs: time
; Output: time in new standard

; Input options:
;  from_ut (also works for FITS formatted dates), from_mjd, from_aft,
;  or from_nustar (MET seconds, default)
; Output options:
;  mjd, fits, aft, nustar (MET seconds, default)
  
  
;; Updated log:
;; 2016/09/06:
;; Removed dependence on cmsysotal and cleaned up code. Added jd option.
;; 2012-08-17
;; Added support for time formats from the AFT, e.g.:
;;               2012:220:00:59:59
;; Brian Grefenstette
; Note: uses date_conv (astrolib) extensively

;  print, time, format = '(d0)'
  if  strcmp(typename(time), 'FLOAT') then time = double(time)
;  print, typename(time)
;  print, time, format = '(d0)'
  mjdref = double(55197.)       ; reference time for NuSTAR. NuSTAR time is seconds
               ; since this epoch

; Get offset for leap seconds:
set = 0
ntime = n_elements(time)

IF keyword_set(from_ut) THEN BEGIN
   ; Converting from UTC string to NuSTAR epoch time:
                                ; Format for time string should be;
;               DD-MON-YEAR HH:MM:SS.SS
;               (eg.  14-JUL-2005 15:25:44.23)
;            OR
;               YYYY-MM-DD HH:MM:SS.SS  (ISO standard)

   ; Convert this to JD:
   jd_out_time = dblarr(ntime)
   FOR i = 0, ntime - 1 DO begin
      jd_out_time[i] = date_conv(time[i], 'J')
   ENDFOR
endif else if keyword_set(from_mjd) THEN BEGIn

   ; Convert to JD from MJD
   jd_out_time = time + 2400000.5

endif else if keyword_set(from_aft) THEN BEGIN
   ; Converting from AFT string to NuSTAR epoch time:
                                ; Format for time string should be;
;                YEAR:DOY:HR:MN:SS.SS
;                2012:220:00:59:59
   jd_out_time = dblarr(ntime)

   FOR i = 0, ntime - 1 DO BEGIN
      temp = float(strsplit(time[i], ':', /extract))
      jd_out_time[i] = date_conv(temp, 'J')
   ENDFOR  
ENDIF else begin

   ; Convert NuSTAR MET seconds to JD:
   mjd_out_time = (double(time) / double(86400.) ) + mjdref
   jd_out_time = double(mjd_out_time)+ double(2400000.5)

endelse






; Now format output:

IF keyword_set(mjd) THEN BEGIN
   output = jd_out_time - double(2400000.5)


   set = 1
ENDIF else IF keyword_set(fits) THEN BEGIN
   
   output = strarr(ntime)
   FOR i = 0, ntime - 1 DO begin
      output[i] =  date_conv(jd_out_time[i], 'F')
   ENDFOR
   set = 1
ENDIF else IF keyword_set(ut) THEN BEGIN
   
   output = strarr(ntime)
;   print, ntime
   FOR i = 0, ntime - 1 DO begin
      tmp =  date_conv(jd_out_time[i], 'F')

      tmp = strsplit(tmp, 'T', /extract)
;      print, tmp
      output[i] = tmp[0] + ' '+tmp[1]
   ENDFOR
ENDIF else if keyword_set(aft) THEN BEGIN
   
   output = strarr(ntime)
;   print, ntime
   FOR i = 0, ntime - 1 DO begin
      tmp =  date_conv(jd_out_time[i], 'V')
      output[i] = string(tmp[0], format = '(i0)')+':'+$
            string(tmp[1], format = '(i0)')+':'+$
            string(tmp[2], format = '(i0)')+':'+$
            string(tmp[3], format = '(i0)')+':'+$
            string(tmp[4], format = '(i0)')
   ENDFOR
   set = 1
ENDIF else if keyword_set(jd) then begin
   output = jd_out_time
endif else begin
   output = (double(jd_out_time) - double(2400000.5) - mjdref) * double(86400.)
endelse


return, output


END

