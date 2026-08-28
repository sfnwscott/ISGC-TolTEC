# load packages
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import QTable
import time
import os

# put in current directory with all the redu__ folders

def pointing_correction(obs_map, full_map, obs_header, full_header):
    # get wcs for obs and full map. clip max and min
    wcs_o = WCS(obs_header,naxis=2)
    wcs_f = WCS(full_header,naxis=2)
    RA, DEC = wcs_f.wcs.crval
    center_coordinate = SkyCoord(ra=RA, dec=DEC, frame='icrs', obstime="J2000", unit=(u.deg,u.deg))
    
    full_cy, full_cx = np.round(center_coordinate.to_pixel(wcs_f))
    full_cy, full_cx = int(full_cy), int(full_cx)
    full_map = full_map[full_cy-180:full_cy+180,full_cx-180:full_cx+180]

    obs_cy, obs_cx = np.round(center_coordinate.to_pixel(wcs_o))
    obs_cy, obs_cx = int(obs_cy), int(obs_cx)
    obs_map = obs_map[obs_cy-80:obs_cy+80,obs_cx-80:obs_cx+80]
    obs_h, obs_w = obs_map.shape

    grid_len_y = 32
    grid_len_x = 32
    x, y = np.meshgrid(np.arange(-grid_len_x//2, grid_len_x//2),np.arange(-grid_len_y//2, grid_len_y//2))
    grid = np.stack((x,y),axis=-1).reshape((grid_len_x*grid_len_y,2))
    cross_corrs, rmses = [], []
    for i,disp in enumerate(grid):
        obs_exp = np.zeros(full_map.shape)
        obs_exp[100+disp[0]:100+obs_h+disp[0],100+disp[1]:100+obs_w+disp[1]] = obs_map
        
        cross_corr_ = np.sum(obs_exp * full_map) / (obs_h*obs_w)
        full_int = full_map[100+disp[0]:100+obs_h+disp[0],100+disp[1]:100+obs_w+disp[1]]
        rmse_ = np.sqrt(np.mean((full_int - obs_map)**2))
        cross_corrs.append(cross_corr_)
        rmses.append(rmse_)

    cross_corrs = np.array(cross_corrs).reshape((len(cross_corrs),1))
    stacked = np.hstack((grid, cross_corrs))
    stat_grid = np.zeros((grid_len_y,grid_len_x))
    stat_grid[stacked[:,0].astype(int)+grid_len_y//2,stacked[:,1].astype(int)+grid_len_x//2] = stacked[:,2]
   
    peak_S = np.max(stat_grid)
    peak_S_y, peak_S_x = np.unravel_index(np.argmax(stat_grid), stat_grid.shape)
    offset_S_y = peak_S_y - grid_len_y//2
    offset_S_x = peak_S_x - grid_len_x//2

    rmses = np.array(rmses).reshape((len(rmses),1))
    stacked = np.hstack((grid, rmses))
    stat_grid = np.zeros((grid_len_y,grid_len_x))
    stat_grid[stacked[:,0].astype(int)+grid_len_y//2,stacked[:,1].astype(int)+grid_len_x//2] = stacked[:,2]
   
    peak_N = np.min(stat_grid)
    peak_N_y, peak_N_x = np.unravel_index(np.argmin(stat_grid), stat_grid.shape)
    offset_N_y = peak_N_y - grid_len_y//2
    offset_N_x = peak_N_x - grid_len_x//2
    return offset_S_y, offset_S_x, peak_S, offset_N_y, offset_N_x, peak_N

time_start = time.time()

directory = os.listdir()
obs_names = [131947,131949,134710,134712,134714,134859,134861]
neigs, obs_name, offsets_y, offsets_x, peaks_S, offsets_N_y, offsets_N_x, peaks_N = [], [], [], [], [], [], [], []
full_path = 'redu02/coadded/raw/toltec_commissioning_a1100_citlali.fits'
full_fits_file = fits.open(full_path)
full_map = full_fits_file[1].data[0][0]
full_header = full_fits_file[1].header
for path in directory:
    if not os.path.isdir(path):
        continue
    else:
        pass
    neig = full_fits_file[0].header['HIERARCH CONFIG.CLEANED.NEIG']
    for obs in obs_names:
        path_obs = path + '/' + str(obs) + '/raw/toltec_commissioning_a1100_science_' + str(obs) + '_citlali.fits'
        obs_fits_file_ = fits.open(path_obs)
        map_dat = obs_fits_file_[1].data[0][0]
        map_header = obs_fits_file_[1].header
        offset_y, offset_x, peak_S, offset_N_y, offset_N_x, peak_N = pointing_correction(map_dat, full_map, map_header, full_header)
        offsets_y.append(offset_y)
        offsets_x.append(offset_x)
        peaks_S.append(float(f'{peak_S:.3e}'))
        offsets_N_y.append(offset_N_y)
        offsets_N_x.append(offset_N_x)
        peaks_N.append(float(f'{peak_N:.3e}'))
        neigs.append(neig)
        obs_name.append(str(obs))
        
pointing_offsets = QTable([neigs, peaks_S, peaks_N, obs_name, offsets_x, offsets_y, offsets_N_x, offsets_N_y], 
                          names=('N-eig', 'Peak S', 'Peak N', 'Obs', 'x Offset', 'y Offset', 'x Offset N', 'y Offset N'))
pointing_offsets.pprint_all()
pointing_offsets.write('signal_noise_table.csv', format='ascii.csv', overwrite=True)

time_end = time.time()
print(f'Time Elapsed: {round(time_end-time_start)}s')



        

