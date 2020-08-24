# Estimate Moon's position on an image with Moon and stars
#
# Peter Monta, August 2020

import numpy as np

from scipy.spatial.transform import Rotation as R
import scipy.signal
import scipy.optimize

from astropy.io import fits
from skyfield.api import load, Topos

import sys
import subprocess

# synthesize a flat-shaded moon.
#   size: size of output image in pixels
#   moon_diameter:  diameter of moon in pixels (should be less than size)
#   sun_angle:  position angle to sun in degrees
#   phase: sun-moon-earth angle in degrees

def synthesize_moon(size,moon_diameter,sun_angle,phase):
  old_seterr = np.seterr(invalid='ignore')
  s = np.array((0,0,1))
  r = R.from_euler('y', phase, degrees=True)
  s = r.apply(s)
  r = R.from_euler('z', sun_angle, degrees=True)
  s = r.apply(s)
  e = size/moon_diameter
  t = np.linspace(-e,e,size)
  grid = np.meshgrid(t,t)
  x,y = grid[0],grid[1]
  z = np.sqrt(1-x*x-y*y)
  dot = s[0]*x + s[1]*y + s[2]*z
  moon = (dot>0).astype('float')
  np.seterr(**old_seterr)
  return moon

# position of Moon from the JPL ephemeris using Skyfield

def skyfield_ra_dec(year,month,day,hour,minute,second,clock_correction):
  t = ts.utc(year,month,day,hour,minute,second)
  t = ts.tt_jd(t.tt + clock_correction/86400)
# fixme: don't hardcode topocentric position; get from fits OBSGEO*
  mylocation = earth + Topos('37.431121 N', '122.183724 W')
  d = mylocation.at(t).observe(moon).apparent()
  ra, dec, distance = d.radec()
  return ra._degrees, dec.degrees

# rough estimate of Moon's position in pixels by correlating with a flat-shaded Moon

def moon_initial(img):
# fixme: need to derive moon diameter, position angle, phase angle from fits wcs and ephemeris
  m = synthesize_moon(256,90,240,70)
  c = scipy.signal.correlate(img,m,mode="same")
  idx = np.argmax(c)
  idx = np.unravel_index(idx, c.shape)
  my,mx = idx[0],idx[1]
  return mx+1,my+1

def edge(img,c,r):
# fixme: don't hardcode position angle
  theta = np.arange(240-90,240+90,1.0)
  a = theta*np.pi/180
  inds = c[1]+r*np.sin(a), c[0]+r*np.cos(a)
  e = scipy.ndimage.map_coordinates(img,inds,order=3,prefilter=False)
  return np.sum(e)

def metric(c,img):
# fixme: use image scale to obtain radius in pixels
  e1 = edge(img,c,46.2)
  e2 = edge(img,c,46.8)
  return -(e1-e2)

# fine-tune the Moon's position with a metric based on contrast near the limb

def moon_fine(img,x,y):
  x0 = (x-1,y-1)
# fixme: use image scale to set spatial scale to ~40 arcsec
  simplex_size = 2.0
  initial_simplex = (
    (x0[0],x0[1]),
    (x0[0]+simplex_size,x0[1]),
    (x0[0],x0[1]+simplex_size)
  )
  res = scipy.optimize.minimize(metric,x0,method='Nelder-Mead',args=img,options={'xatol':0.0001,'initial_simplex':initial_simplex})
  return res.x[0]+1, res.x[1]+1

# call wcstools's xy2sky to translate pixels to equatorial coordinates; astropy doesn't handle SCAMP's TPV

def xy2sky(filename,x,y):
  cmd = ['xy2sky','-n','8','-d',filename,'%f'%x,'%f'%y]
#  print('xy2sky args: ',cmd)
  output = subprocess.check_output(cmd)
#  print('xy2sky output: ',output)
  x = output.split()
  ra = float(x[0])
  dec = float(x[1])
  return ra, dec

def fits_time(filename):
  pass
# read FITS header, parse time keyword(s)

#
# moon-position.py image.fits
#

filename = sys.argv[1]
year,month,day,hour,minute,second = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]),float(sys.argv[7])
#year,month,day,hour,minute,second = fits_time(filename)

planets = load('de421.bsp')
earth,moon = planets['earth'],planets['moon']
ts = load.timescale()

# fixme: get clock correction from command-line argument? FITS header? or correct FITS time as part of image calibration
clock_correction = -95.49
ra_true, dec_true = skyfield_ra_dec(year,month,day,hour,minute,second,clock_correction)

img = fits.getdata(filename).astype('float')

x,y = moon_initial(img)
x,y = moon_fine(img,x,y)
ra_est, dec_est = xy2sky(filename,x,y)

print(filename, ra_true, dec_true, ra_est, dec_est, ra_est-ra_true, dec_est-dec_true)
