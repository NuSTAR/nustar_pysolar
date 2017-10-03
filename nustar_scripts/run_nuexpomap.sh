#!/bin/bash

# usage: run_nuexpomap.sh OBSID ENERGY
# e.g.:
# run_products.sh 60001043006 A source.reg bkg.reg
# Dumps data products into a new ./OBSID/source directory


# Set up your local NuSTAR science environment here:
if [ -z "$NUSTARSETUP" ]; then 
    echo "Need to set the NUSTARSETUP environment variable!"
    exit
fi
source $NUSTARSETUP

OBSID=$1
ENERGY=$2

# Modify this to be where you want the maps to land...
OUTDIR=$OBSID/event_cl/expo_maps
if [ ! -d $OUTDIR ] ; then
    mkdir $OUTDIR
fi

for MOD in A B
do
    INSTRUMENT=FPM${MOD}
    DATPATH=${OBSID}/event_cl
    echo $OBSID $MOD $ENERGY

    STEM=nu`basename $OBSID`
    infile=${DATPATH}/${STEM}${MOD}01_cl.evt
    outname=${STEM}${MOD}

    attfile=${DATPATH}/${STEM}_att.fits
    mastaspectfile=${DATPATH}/${STEM}_mast.fits
    det1reffile=${DATPATH}/${STEM}${MOD}_det1.fits

    expomapfile=${OUTDIR}/${outname}01_ex_$ENERGY.img
    expolog=${OUTDIR}/${outname}01_ex_${ENERGY}.log

    echo PID $$ > $expolog

    nuexpomap infile=$infile attfile=$attfile \
        mastaspectfile=$mastaspectfile \
        det1reffile=$det1reffile \
        vignflag=yes \
        energy=$ENERGY \
        pixbin=5 \
        clobber=yes \
        expomapfile=$expomapfile >> $expolog 2>&1
done

