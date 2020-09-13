#!/bin/bash

# Process a raw camera image into a FITS image with WCS.  Uses a master distortion WCS (to avoid overfitting
# when only a few stars are available) and a metadata JSON file containing time and topocentric location.

RAW_FILE=$1
METADATA_FILE=$2
DISTORTION_FILE=$3
FITS_FILE=$4

# process the RAW camera image into a FITS file.  No deBayering.
# fixme: bad pixels, flat fielding

raw2fits.py ${RAW_FILE} ${METADATA_FILE} temp1.fits

# extract green channel, turning 45 degrees (quincunx --> rectangular)

decimate45.py temp1.fits temp2.fits

# Astromatic's source extractor

sex temp2.fits -c sex.txt

# apply a SIP distortion polynomial to the detected sources

predistort.py test.cat predistort.cat ${DISTORTION_FILE}

# transform the FITS_LDAC source catalog into a "xylist" for use with astrometry.net tools

ldac2xylist.py predistort.cat predistort.xylist

# initial astrometry.net solve-field

#fixme: don't hardcode the image size and scale?
solve-field --no-plots --overwrite --x-column XWIN_IMAGE --y-column YWIN_IMAGE --sort-column FLUX_AUTO --width 3283 --height 3283 --scale-low 17 --scale-high 22 --scale-units arcsecperpix --no-tweak --crpix-center predistort.xylist

# save the initial WCS

mv predistort.xylist.wcs init.wcs
cp predistort.xylist predistort2.xylist

# second call to solve-field, to refine the WCS using maximum number of matched sources

#fixme: don't hardcode the image size and scale?
solve-field --verify init.wcs --no-plots --overwrite --x-column XWIN_IMAGE --y-column YWIN_IMAGE --sort-column FLUX_AUTO --width 3283 --height 3283 --scale-low 17 --scale-high 22 --scale-units arcsecperpix --tweak-order 1 --crpix-center predistort2.xylist

# extract the WCS header

fitshdr predistort2.xylist.wcs >temp3.ahead

# make a symbolic link to the predistorted catalog

rm -f temp3.cat
ln -s predistort.cat temp3.cat

# run SCAMP on the catalog and astrometry.net's WCS, making a new refined WCS.  Linear distortion terms only
# (higher-order distortion has already been removed fro the source catalog).

scamp temp3.cat -c scamp.txt

# prepare output file

cp temp2.fits ${FITS_FILE}
cp temp3.head `basename ${FITS_FILE} .fits`.head

# copy the WCS header into the output FITS file

missfits ${FITS_FILE}
