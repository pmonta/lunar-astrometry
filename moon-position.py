# Estimate Moon's position on an image with Moon and stars
#
# Peter Monta, August 2020

import numpy as np
import numpy.linalg

from scipy.spatial.transform import Rotation as R
import scipy.signal
import scipy.optimize

from astropy.io import fits
import astropy.time

from skyfield.api import Loader, Topos
from skyfield.constants import AU_M
from skyfield.earthlib import reverse_terra
from skyfield.timelib import Time
from skyfield.trigonometry import position_angle_of
from skyfield.units import Angle

import sip

import sys
import subprocess

load = Loader('~/.skyfield-data', verbose=False)
planets = load('de421.bsp')
earth,moon,sun = planets['earth'],planets['moon'],planets['sun']
ts = load.timescale(builtin=False)

# synthesize a flat-shaded moon.
#   size: size of output image in pixels
#   moon_diameter:  diameter of moon in pixels (should be less than size)
#   sun_angle:  position angle to sun in degrees
#   phase: sun-moon-earth angle in degrees

def synthesize_moon(size,moon_diameter_pixels,sun_angle,moon_phase):
  old_seterr = np.seterr(invalid='ignore')
  s = np.array((0,0,1))
  r = R.from_euler('y', moon_phase, degrees=True)
  s = r.apply(s)
  r = R.from_euler('z', sun_angle, degrees=True)
  s = r.apply(s)
  e = size/moon_diameter_pixels
  t = np.linspace(-e,e,size)
  grid = np.meshgrid(t,t)
  x,y = grid[0],grid[1]
  z = np.sqrt(1-x*x-y*y)
  dot = s[0]*x + s[1]*y + s[2]*z
  moon = (dot>0).astype('float')
  np.seterr(**old_seterr)
  return moon

# Skyfield time and location

def skyfield_time(year,month,day,hour,minute,second):
  return ts.utc(year,month,day,hour,minute,second)

def skyfield_location_xyz(x,y,z):
  lat,lon,e = reverse_terra((x/AU_M,y/AU_M,z/AU_M),gast=0)
  return Topos(latitude_degrees=np.rad2deg(lat), longitude_degrees=np.rad2deg(lon), elevation_m=e)

# radius of Moon in degrees

def skyfield_moon_radius(t,location):
# fixme: incorporate flattening, refraction?
  MOON_RADIUS_KM = 1737.4
  d = location.at(t).observe(moon).apparent()
  ra, dec, distance = d.radec()
  theta = np.arctan(MOON_RADIUS_KM/distance.km)
  return np.rad2deg(theta)

# position of moon from ephemeris

def skyfield_moon_ra_dec(t,location):
  d = location.at(t).observe(moon).apparent()
  ra, dec, distance = d.radec()
  return ra._degrees, dec.degrees

# position of sun from ephemeris

def skyfield_sun_ra_dec(t,location):
  d = location.at(t).observe(sun).apparent()
  ra, dec, distance = d.radec()
  return ra._degrees, dec.degrees

# moon phase (sun-moon-earth)

def skyfield_moon_phase(t,location):
  d1 = moon.at(t).observe(location)
  d2 = moon.at(t).observe(sun)
  return d1.separation_from(d2).degrees

# angle from center of image to Sun

def skyfield_sun_position_angle(filename,t,location):
  hdulist = fits.open(filename)
  hdr = hdulist[0].header
  cx = hdr['NAXIS1']/2  # fixme: use IMAGE{H,W}?
  cy = hdr['NAXIS2']/2
  hdulist.close()
  center_ra,center_dec = xy2sky(filename,cx,cy)
  unit_ra,unit_dec = xy2sky(filename,cx+1,cy)    # fixme: use a true vector, not a finite difference
  sun_ra,sun_dec = skyfield_sun_ra_dec(t,location)
  return pa(center_ra,center_dec,sun_ra,sun_dec) - pa(center_ra,center_dec,unit_ra,unit_dec)

def pa(c_ra,c_dec,a_ra,a_dec):
  c_ra = Angle(hours=c_ra/15)
  c_dec = Angle(degrees=c_dec)
  a_ra = Angle(hours=a_ra/15)
  a_dec = Angle(degrees=a_dec)
  return position_angle_of((c_dec,c_ra),(a_dec,a_ra)).degrees

# rough estimate of Moon's position in pixels by correlating with a flat-shaded Moon

def moon_initial(img,t,location,plate_scale,moon_radius,sun_position_angle):
  moon_diameter_pixels = 2*moon_radius*3600/plate_scale
  size = int(moon_diameter_pixels*1.2)
  moon_phase = skyfield_moon_phase(t,location)
  m = synthesize_moon(size,moon_diameter_pixels,sun_position_angle,moon_phase)
  c = scipy.signal.correlate(img,m,mode="same")
  idx = np.argmax(c)
  idx = np.unravel_index(idx, c.shape)
  my,mx = idx[0],idx[1]
  return mx+1,my+1

# intensity along lunar limb

def edge(img,c,r):
# fixme: use only 180 degrees of limb
  theta = np.arange(0,360,1.0)
  a = np.deg2rad(theta)
  inds = c[1]+r*np.sin(a), c[0]+r*np.cos(a)
  e = scipy.ndimage.map_coordinates(img,inds,order=3,prefilter=False)
  return np.sum(e)

# optimizer metric: difference between edge(r) and edge(r+dr)

def metric(c,img,moon_radius_pixels):
  p = (c[0],c[1])
  r = moon_radius_pixels
  dr = 0.6
  e1 = edge(img,p,r-dr/2)
  e2 = edge(img,p,r+dr/2)
  return -(e1-e2)

# fine-tune the Moon's position with a metric based on contrast near the limb

def moon_fine(img,x,y,plate_scale,moon_radius):
  moon_radius_pixels = moon_radius*3600/plate_scale
  moon_radius_pixels += 0.4
  x0 = (x-1,y-1)
  simplex_size = 40.0/plate_scale
  initial_simplex = (
    (x0[0],x0[1]),
    (x0[0]+simplex_size,x0[1]),
    (x0[0],x0[1]+simplex_size)
  )
  res = scipy.optimize.minimize(metric,x0,method='Nelder-Mead',args=(img,moon_radius_pixels),options={'xatol':0.0001,'initial_simplex':initial_simplex})
  return res.x[0]+1, res.x[1]+1

# call wcstools's xy2sky to translate pixels to equatorial coordinates; astropy doesn't handle SCAMP's TPV

def xy2sky(filename,x,y):
  cmd = ['xy2sky','-n','12','-d',filename,'%f'%x,'%f'%y]
  output = subprocess.check_output(cmd)
  x = output.split()
  ra = float(x[0])
  dec = float(x[1])
  return ra, dec

# observation time from FITS header

def fits_time(filename):
  hdulist = fits.open(filename)
  hdr = hdulist[0].header
  if 'TIMESYS' in hdr and hdr['TIMESYS']!='UTC':
    raise 'TIMESYS must be UTC'
  iso = hdr['DATE-OBS']
  hdulist.close()
  t = astropy.time.Time(iso, format='isot', scale='utc')
  d = t.datetime
  year,month,day = d.year,d.month,d.day
  hour,minute,second = d.hour,d.minute,d.second
  microsecond = d.microsecond
  return year, month, day, hour, minute, second+microsecond/1000000.0

# location (ITRF) from FITS header

def fits_xyz(filename):
  hdulist = fits.open(filename)
  hdr = hdulist[0].header
  x = hdr['OBSGEO-X']
  y = hdr['OBSGEO-Y']
  z = hdr['OBSGEO-Z']
  hdulist.close()
  return x,y,z

# calculate average plate scale in arcseconds per pixel

def fits_plate_scale(filename):
  hdulist = fits.open(filename)
  hdr = hdulist[0].header
  cd11 = hdr['CD1_1']
  cd12 = hdr['CD1_2']
  cd21 = hdr['CD2_1']
  cd22 = hdr['CD2_2']
  hdulist.close()
  cd = np.array([[cd11,cd12],[cd21,cd22]])
  scale = np.sqrt(np.abs(np.linalg.det(cd)))
  return scale*3600

def project(ra1,dec1,ra2,dec2):
  ra1,dec1 = np.deg2rad(ra1),np.deg2rad(dec1)
  ra2,dec2 = np.deg2rad(ra2),np.deg2rad(dec2)
  s = np.array((np.cos(ra2)*np.cos(dec2), np.sin(ra2)*np.cos(dec2), np.sin(dec2)))
  r = R.from_euler('z', -ra1)
  s = r.apply(s)
  r = R.from_euler('y', dec1)
  s = r.apply(s)
  xi = np.rad2deg(np.arcsin(s[1]))
  eta = np.rad2deg(np.arcsin(s[2]))
  return xi, eta

#
# moon-position.py image.fits sip.fits
#

image_filename = sys.argv[1]
sip_filename = sys.argv[2]

year,month,day,hour,minute,second = fits_time(image_filename)
geo_x,geo_y,geo_z = fits_xyz(image_filename)

t = skyfield_time(year,month,day,hour,minute,second)
location = earth + skyfield_location_xyz(geo_x,geo_y,geo_z)
moon_radius = skyfield_moon_radius(t,location)
ra_true, dec_true = skyfield_moon_ra_dec(t,location)
sun_position_angle = skyfield_sun_position_angle(image_filename,t,location)
plate_scale = fits_plate_scale(image_filename)

img = fits.getdata(image_filename).astype('float')

x,y = moon_initial(img,t,location,plate_scale,moon_radius,sun_position_angle)
x,y = moon_fine(img,x,y,plate_scale,moon_radius)

xpd,ypd = sip.sip(sip_filename,x,y)
ra_est, dec_est = xy2sky(image_filename,xpd,ypd)

xi, eta = project(ra_true,dec_true,ra_est,dec_est)

print("%s %f %f %f %f %.3f %.3f" % (image_filename, ra_true, dec_true, ra_est, dec_est, 3600*xi, 3600*eta))
