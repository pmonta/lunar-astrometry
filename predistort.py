#!/usr/bin/python

from astropy.io import fits
import sys

#fixme: import sip

def read_sip(filename):
  hdu_list = fits.open(filename)
  return hdu_list[0].header

def sip_eval_a(sip,u,v):
  a00 = sip['A_0_0']
  a01 = sip['A_0_1']
  a02 = sip['A_0_2']
  a03 = sip['A_0_3']
  a10 = sip['A_1_0']
  a11 = sip['A_1_1']
  a12 = sip['A_1_2']
  a20 = sip['A_2_0']
  a21 = sip['A_2_1']
  a30 = sip['A_3_0']
  x = a00 + a01*v + a02*v*v + a03*v*v*v + a10*u + a11*u*v + a12*u*v*v + a20*u*u + a21*u*u*v + a30*u*u*u
  return x

def sip_eval_b(sip,u,v):
  b00 = sip['B_0_0']
  b01 = sip['B_0_1']
  b02 = sip['B_0_2']
  b03 = sip['B_0_3']
  b10 = sip['B_1_0']
  b11 = sip['B_1_1']
  b12 = sip['B_1_2']
  b20 = sip['B_2_0']
  b21 = sip['B_2_1']
  b30 = sip['B_3_0']
  x = b00 + b01*v + b02*v*v + b03*v*v*v + b10*u + b11*u*v + b12*u*v*v + b20*u*u + b21*u*u*v + b30*u*u*u
  return x

def predistort(x,y,cx,cy,sip):
  ex = sip_eval_a(sip,x-cx,y-cy)
  ey = sip_eval_b(sip,x-cx,y-cy)
  return x+ex, y+ey

#
# predistort.py input.fits output.fits sip.fits
#

input_filename = sys.argv[1]
output_filename = sys.argv[2]
sip_filename = sys.argv[3]

sip = read_sip(sip_filename)

hdu_list = fits.open(input_filename)
data = hdu_list[2].data

#fixme
#cx = ['CRPIX1']
#cy = ['CRPIX2']
cx = 1642
cy = 1642
data['XWIN_IMAGE'],data['YWIN_IMAGE'] = predistort(data['XWIN_IMAGE'],data['YWIN_IMAGE'],cx,cy,sip)

hdu_list.writeto(output_filename,overwrite=True)
