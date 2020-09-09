from astropy.io import fits
from astropy.time import Time, TimeDelta

import exifread
import datetime
import subprocess
import json
import sys

def dcraw(image_filename, pgm_filename):
  cmd = ['dcraw','-t','0','-E','-4','-j','-c',image_filename]
  with open(pgm_filename,"w") as fp_out:
    subprocess.run(cmd, stdout=fp_out)

def pgm_to_fits(pgm_filename, fits_filename):
# fixme: do in pure python
  cmd = ['pamflip','-topbottom',pgm_filename]
  with open('temp2.pgm',"w") as fp_out:
    subprocess.run(cmd, stdout=fp_out)
  cmd = ['pamtofits','temp2.pgm']
  with open(fits_filename,"w") as fp_out:
    subprocess.run(cmd, stdout=fp_out)

def datetime_iso(exif_datetime,exif_subsec,clock_correction):
  x = exif_datetime.values.split()
  d = x[0].split(':')
  year, month, day = int(d[0]), int(d[1]), int(d[2])
  e = x[1].split(':')
  hour, minute, second = int(e[0]), int(e[1]), int(e[2])
  t = Time(datetime.datetime(year, month, day, hour, minute, second), scale='utc')
  t += TimeDelta(int(exif_subsec.values)/100.0, format='sec')
  t += TimeDelta(clock_correction, format='sec')
  return t.fits

#
# raw2fits.py image.cr2 metadata.json image.fits
#

image_filename = sys.argv[1]
json_filename = sys.argv[2]
fits_filename = sys.argv[3]
pgm_filename = 'temp.pgm'

dcraw(image_filename, pgm_filename)

pgm_to_fits(pgm_filename, 'temp.fits')

hdu_list = fits.open('temp.fits')

with open(json_filename) as fp:
  m = json.load(fp)

with open(image_filename,"rb") as exif_fp:
  exif_tags = exifread.process_file(exif_fp, details=False)

hdr = hdu_list[0].header

hdr['DATE-OBS'] = datetime_iso(exif_tags['EXIF DateTimeOriginal'],exif_tags['EXIF SubSecTimeOriginal'],m['clock-correction'])
hdr['OBSGEO-X'] = m['location-ITRF']['x']
hdr['OBSGEO-Y'] = m['location-ITRF']['y']
hdr['OBSGEO-Z'] = m['location-ITRF']['z']

# fixme: additional header info
#hdr['DATE']    ==date-obs
#hdr['TIMESYS'] fixed, 'UTC'
#hdr['XPOSURE'] exif
#hdr['GAIN']     json
#hdr['INSTRUME'] json
#    comment: Camera Model Name       (exif)
#    comment: camera serial number    (exif)
#    comment: lens                    (json)
#    comment: ISO                     (exif)
#    comment: sensor temperature      (exif)

hdu_list.writeto(fits_filename,overwrite=True)
