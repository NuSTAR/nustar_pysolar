import numpy as np
import astropy.time


def bad_pix(evtdata, fpm='A'):
    """Do some basic filtering on known bad pixels.
    
    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"FPMA" | "FPMB"}
        Which FPM you're filtering on. Assumes A if not set.


    Returns
    -------

    goodinds: iterable
        Index of evtdata that passes the filtering.
    """
    

    
    # Hot pixel filters
    
    # FPMA or FPMB
    
    if fpm.find('B') == -1 :
        pix_filter = np.invert( ( (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 16) & (evtdata['RAWY'] == 5) |
                                (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 22) |
                                (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 27) & (evtdata['RAWY'] == 6) |
                                (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 27) & (evtdata['RAWY'] == 21) |
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 22) & (evtdata['RAWY'] == 1) |
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 15) & (evtdata['RAWY'] == 3) |
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 5) & (evtdata['RAWY'] == 5) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 22) & (evtdata['RAWY'] == 7) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 16) & (evtdata['RAWY'] == 11) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 18) & (evtdata['RAWY'] == 3) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 4) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 25) & (evtdata['RAWY'] == 5) ) )
    else:
        pix_filter = np.invert( ( (evtdata['DET_ID'] == 0) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 24)) )


    inds = (pix_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds
    
def by_energy(evtdata, energy_low=2.5, energy_high=10.):
    """ Apply energy filtering to the data.
    
    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.
        
    energy_low: float
        Low-side energy bound for the map you want to produce (in keV).
        Defaults to 2.5 keV.

    energy_high: float
        High-side energy bound for the map you want to produce (in keV).
        Defaults to 10 keV.
    """     
    pilow = (energy_low - 1.6) / 0.04
    pihigh = (energy_high - 1.6) / 0.04
    pi_filter = ( ( evtdata['PI']>pilow ) &  ( evtdata['PI']<pihigh))
    inds = (pi_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds

def by_det(evtdata, det_id=0):
    """ Apply Det filtering to the data.
    
    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.
        
    det_id: int 
        Values of 0,1,2 or 3
        Defaults to 0


    """     
    det_filter = (evtdata['DET_ID'] == det_id)
    inds = (det_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds


def by_xy(evtdata, hdr, xy_range):
    """ Apply position filtering to the data.
    
    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.
        
    hdr: FITS header structure
        This should be an hdu.header structure from a NuSTAR FITS file
        
    xy_range: float list 4 elements
        Min, max x and min, max y in HPS arcseconds
    """     
    
#   Convert the xy HPS arcseconds limits into NuSTAR pixel locations
    for field in hdr.keys():
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
                xval = field[5:8]
            if hdr[field] == 'Y':
                yval = field[5:8]
#   x,y should be the same but just incase
    npixx=hdr['TLMAX'+xval]
    npixy=hdr['TLMAX'+yval]
    pixsizex=np.abs(hdr['TCDLT'+xval])
    pixsizey=np.abs(hdr['TCDLT'+yval])
    
    min_x=(xy_range[0]/pixsizex)+npixx*0.5
    max_x=(xy_range[1]/pixsizex)+npixx*0.5
    
    min_y=(xy_range[2]/pixsizey)+npixy*0.5
    max_y=(xy_range[3]/pixsizey)+npixy*0.5

    xy_filter = ( ( evtdata['X']>min_x ) &  ( evtdata['X']<max_x) & \
                ( evtdata['Y']>min_y ) &  ( evtdata['Y']<max_y))
    inds = (xy_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds

def by_time(evtdata, hdr, time_range):
    """ Apply time filtering to the data.
    
    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.
        
    hdr: FITS header structure
        This should be an hdu.header structure from a NuSTAR FITS file
        
    time_range: astropy.time format 2 elements
        Min, max time range to consider, in astropy time format
    """     
    
#   Convert the time_range into seconds from NuSTAR ref time - units of evtdata['time']
    ns_ref=astropy.time.Time(hdr['MJDREFI'],format='mjd')
    
    min_t=(time_range[0].unix-ns_ref.unix)
    max_t=(time_range[1].unix-ns_ref.unix)
#     print(min_t,max_t)
    time_filter = ( ( evtdata['time']>min_t ) &  ( evtdata['time']<max_t))
    inds = (time_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds
    
def gradezero(evtdata):
    """ Only accept counts with GRADE==0.
        
    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.
        
    Returns
    -------

    goodinds: iterable
        Index of evtdata that passes the filtering.
    """

    # Grade filter
    
    grade_filter = ( evtdata['GRADE'] == 0)
    inds = (grade_filter).nonzero()
    goodinds = inds[0]
    
    return goodinds

def event_filter(evtdata, fpm='FPMA',
    energy_low=2.5, energy_high=10,hdr=0,xy_range=0,time_range=0,
                 dets_id=[],no_grade_filter=False,no_bad_pix_filter=False):
    """ All in one filter module. By default applies an energy cut, 
        selects only events with grade == 0, and removes known hot pixel.
        
        Note that this module returns a cleaned eventlist rather than
        the indices to the cleaned events.

    Parameters
    ----------
    evtdata: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"FPMA" | "FPMB"}
        Which FPM you're filtering on. Defaults to FPMA.
    
    energy_low: float
        Low-side energy bound for the map you want to produce (in keV).
        Defaults to 2.5 keV.

    energy_high: float
        High-side energy bound for the map you want to produce (in keV).
        Defaults to 10 keV.
        
    hdr: FITS header structure
        This should be an hdu.header structure from a NuSTAR FITS file
        Only needed when doing filtering by time or xy
        If not supplied no tim or xy filtering, even if they are specified
        
    xy_range: float list 4 elements
        Min, max x and min, max y in HPS arcseconds
        Default not used, needs hdr to be supplied
        
    time_range: astropy.time format 2 elements
        Min, max time range to consider, in astropy time format
        Default not used, needs hdr to be supplied 
        
    dets_id: int list up to 4 elements
        Which dets (0,1,2 and/or 3) to include
        Defaults to using all dets
        
    no_bad_pix_filter: bool True or False
        Don't filter the bad pixels
        Default is False (so does filter bad pixels)  
        
    no_grade_filter: bool True or False
        Return all grades, not just grade == 0
        Default is False (so does return grade ==0)

    Returns
    -------

    cleanevt: FITS data class.
        This is the subset of evtdata that pass the data selection cuts.
    """
#   Option to not filter out the bad pixels
    if no_bad_pix_filter == False:
        goodinds = bad_pix(evtdata, fpm=fpm)
        evt_badfilter = evtdata[goodinds]
    else:
        evt_badfilter = evtdata
    
#   Filter by energy
    goodinds = by_energy(evt_badfilter,
                         energy_low=energy_low, energy_high = energy_high)
    evt_energy = evt_badfilter[goodinds]
    
#   Option to not filter the grades
    if no_grade_filter == False:
        goodinds = gradezero(evt_energy)
        cleanevt = evt_energy[goodinds]
    else:
        cleanevt = evt_energy 
    
#   Filter for specific dets
    if (len(dets_id) > 0):
        goodinds=np.array([],dtype='int')
        for dd in dets_id:
            gi_temp=by_det(cleanevt,det_id=dd)
            goodinds=np.append(goodinds,gi_temp)
        cleanevt=cleanevt[goodinds]
    
#   Filters for time and postion 
    if hdr != 0 :
#       Only do xy filtering if specified
        if (xy_range != 0):
            goodinds=by_xy(cleanevt,hdr,xy_range)
            cleanevt=cleanevt[goodinds]
#       Only do time filtering if specified
        if (time_range != 0):
            goodinds=by_time(cleanevt,hdr,time_range)
            cleanevt=cleanevt[goodinds]
            
    return cleanevt


