def download_occultation_times(outdir='./'):
    """
    
    Pull the orbital information from the Berkeley website
    and then returns the location of the file.
     
    Parameters
    ----------
    
    outdir: Optional desired output location. Defaults to the working directory.
    
    Returns
    ----------
    
    Returns the filename that you've downloaded.
    
    Notes
    ---------
    
    Will not work if year_day is in the future. This is the time the file was
    generated, not the date you want to observe. The occultation windows should
    extend roughly a month into the future, though they may be more uncertain
    the further ahead you go.
    
    """
    import os
    from datetime import date, timedelta


    # Make sure you've got a trailing slash...
    if not(outdir.endswith('/')):
        outdir+'/'
    
    # Make sure the directory exists and create one if not.
    directory = os.path.dirname(outdir)
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    # Get yesterday's file...
    today = date.today() - timedelta(1) 
    year = str(today.timetuple().tm_year)
    day='{0:03d}'.format(today.timetuple().tm_yday)
    year_doy=year+'_'+day

    
    myname='nustar_pysolar.planning'
    import wget
    
    url='https://ops-srvr.ssl.berkeley.edu/ground_systems/products//NUSTAR/'
    url+=year_doy+'/'
    fname='NUSTAR.'+year_doy+'.SHADOW_ANALYSIS.txt' 
    url+=fname

    # Check to see if the file exists:
    import os.path
    if not(os.path.isfile(outdir+fname)):
        wget.download(url, out=outdir+fname)
    
    return outdir+fname

def download_tle(outdir='./'):
    """Download the NuSTAR TLE archive.
    
    Parameters
    ----------
    
    outdir: Optional desired output location. Defaults to the working directory.
    
    Returns
    ----------
    
    Returns the filename that you've downloaded.
    
    Notes
    ---------
    
    """
    import os
    import wget


    # Make sure you've got a trailing slash...
    if not(outdir.endswith('/')):
        outdir+'/'
    
    # Make sure the directory exists and create one if not.
    directory = os.path.dirname(outdir)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    myname='nustar_pysolar.io.download_tle'
    
    url='http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/NuSTAR.tle'

    # Check to see if the file exists:
    fname = 'NuSTAR.tle'
    outfile = outdir+'/'+fname
    if (os.path.isfile(outfile)):
        os.remove(outfile)
    
    
    wget.download(url, out=outfile)
    
    
    return outfile

def read_tle_file(tlefile, **kwargs):
    """Read in the TLE file. 
    
    Returns the times for each line element and the TLE time
    
    """
    times = []
    line1 = []
    line2 = []
    
    from datetime import datetime
    # Catch if the file can't be opened:
    try:
        f = open(tlefile, 'r')
    except FileNotFoundError:
        print("Unable to open: "+tlefile)
    
    ln=0
    for line in f:
#        print(line)
        if (ln == 0):
            year= int(line[18:20])
            day = int(line[20:23])
            
            times.extend([datetime.strptime("{}:{}".format(year, day), "%y:%j")])
            line1.extend([line.strip()])
            ln=1
        else:
            ln=0
            line2.extend([line.strip()])
    f.close()
    return times, line1, line2

def get_epoch_tle(epoch, tlefile):
    """Find the TLE that is closest to the epoch you want to search.
    
    epoch is a datetime object, tlefile is the file you want to search through.
    
    """

    times, line1, line2 = read_tle_file(tlefile)
    from astropy.time import Time
    
    # Allow astropy Time objects
    if type(epoch) is Time:
        epoch = epoch.datetime
        
    mindt = 100.
    min_ind = 0
    for ind, t in enumerate(times):
        dt = abs((epoch -t).days)
        if dt < mindt:
            min_ind = ind
            mindt = dt

    good_line1 = line1[min_ind]
    good_line2 = line2[min_ind]

    return mindt, good_line1, good_line2
