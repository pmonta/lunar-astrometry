from astropy.io import fits
import sys

#
# ldac2xylist.py input.fits output.fits
#

input_filename = sys.argv[1]
output_filename = sys.argv[2]

hdu_list = fits.open(input_filename)
h = hdu_list[2]

p = fits.PrimaryHDU()
hdu_list_out = fits.HDUList([p,h])

hdu_list_out.writeto(output_filename,overwrite=True)
