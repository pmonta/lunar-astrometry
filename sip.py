from astropy.io import fits
import sys

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

def sip_apply(x,y,cx,cy,sip):
  ex = sip_eval_a(sip,x-cx,y-cy)
  ey = sip_eval_b(sip,x-cx,y-cy)
  return x+ex, y+ey

def sip(filename,x,y):
  sip_header = read_sip(filename)
  cx = sip_header['CRPIX1']
  cy = sip_header['CRPIX2']
  x,y = sip_apply(x,y,cx,cy,sip_header)
  return x,y
