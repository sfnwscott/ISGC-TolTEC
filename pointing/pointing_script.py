# load packages
import numpy as np
import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.modeling import models, fitting
from astroquery.simbad import Simbad
from netCDF4 import Dataset

def locate_pointsource(map,x0_f,y0_f):
    x0_f, y0_f = int(np.round(x0_f)), int(np.round(y0_f))
    map = map[y0_f-20:y0_f+20,x0_f-20:x0_f+20]
    amp = np.max(map)
    y0, x0 = np.unravel_index(map.argmax(), map.shape)
    gauss = models.Gaussian2D(amp, x0, y0, 5, 5, 0)
    fitter = fitting.LevMarLSQFitter(calc_uncertainties=False)
    yi, xi = np.indices(map.shape)
    fit_func = fitter(gauss, xi, yi, map)
    x = fit_func.x_mean + x0_f - 20
    y = fit_func.y_mean + y0_f - 20
    return (x, y)

def rot_a(angle):
    m = np.eye(2)
    m[0,0] = m[1,1] = np.cos(angle)
    m[0,1] = -np.sin(angle)
    m[1,0] = +np.sin(angle)
    return m

path = './pointing_maps/cleaned_maps/'
obs_names = [131947,131949,134710,134712,134714,134859,134861]
simbad_results = Simbad.query_object('VLA 1623')
RA, DEC = simbad_results.to_pandas()[['RA','DEC']].values[0]
real_coordinate = SkyCoord(ra=RA, dec=DEC, frame='icrs', obstime="J2000", unit=(u.hourangle,u.deg))
nc_files = os.listdir('./pointing/nc_files')

for obs in obs_names:
    fits_file_ = fits.open(path + f'toltec_commissioning_a1100_science_{str(obs)}_citlali.fits')
    map_dat = fits_file_[6].data[0][0]
    map_header = fits_file_[6].header
    nc_file_path = [file for file in nc_files if str(obs) in file][0]
    nc_file = Dataset('./pointing/nc_files/' + nc_file_path)
    Tel_PA = nc_file.variables['Data.TelescopeBackend.ActParAng'][:].data # in radians
    Tel_PA = np.median(Tel_PA)
    coord_pix = WCS(map_header,naxis=2).world_to_pixel(real_coordinate)
    mask = map_dat > 100
    cut_map = map_dat * mask
    map_x, map_y = locate_pointsource(cut_map, coord_pix[0], coord_pix[1])
    map_coordinate = WCS(map_header,naxis=2).pixel_to_world(map_x, map_y).icrs # in ICRS
    delta_ra = map_coordinate.ra.arcsec - real_coordinate.ra.arcsec
    delta_dec = map_coordinate.dec.arcsec - real_coordinate.dec.arcsec
    delta_az, delta_alt = rot_a(Tel_PA) @ np.array([delta_ra, delta_dec])
    print(delta_dec, delta_ra)
    
    
    

 
