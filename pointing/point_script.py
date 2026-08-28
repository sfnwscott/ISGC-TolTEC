import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import QTable

def stack_fg(data):
    viewdict = {'signal':1,'weight':2,'kernel':3,'coverage':4,'coverage_bool':5,'sig2noise':6}
    views = []
    for view in ['signal']:
        inds = [viewdict[view] + _ for _ in [0,6,12,18,24]]
        hdu_dict = {}
        for i in inds:
            data_ = data[i].data[0][0]
            hdu_dict[i] = data_
        view_data = sum(np.array(list(hdu_dict.values())))/len(inds)
        views.append(view_data)
    return views

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

    grid_len_y = 36
    grid_len_x = 36
    x, y = np.meshgrid(np.arange(-grid_len_x//2, grid_len_x//2),np.arange(-grid_len_y//2, grid_len_y//2))
    grid = np.stack((x,y),axis=-1).reshape((grid_len_x*grid_len_y,2))
    cross_corrs = []
    for i,disp in enumerate(grid):
        obs_exp = np.zeros(full_map.shape)
        obs_exp[100+disp[0]:100+obs_h+disp[0],100+disp[1]:100+obs_w+disp[1]] = obs_map
        
        # power of 4 to emphasize bright points as much as possible
        cross_corr_ = np.sum(obs_exp**4 * full_map**4) / (obs_h*obs_w)
        cross_corrs.append(cross_corr_)

    cross_corrs = np.array(cross_corrs).reshape((len(cross_corrs),1))
    stacked = np.hstack((grid, cross_corrs))
    stat_grid = np.zeros((grid_len_y,grid_len_x))
    stat_grid[stacked[:,0].astype(int)+grid_len_y//2,stacked[:,1].astype(int)+grid_len_x//2] = stacked[:,2]
   
    peak_S = np.max(stat_grid)
    peak_S_y, peak_S_x = np.unravel_index(np.argmax(stat_grid), stat_grid.shape)
    offset_S_y = peak_S_y - grid_len_y//2
    offset_S_x = peak_S_x - grid_len_x//2
    # sorted_indices = np.argsort(stacked[:, 2])[::-1]
    # print(stacked[sorted_indices][:5])
    return offset_S_y, offset_S_x, peak_S

# read in fits maps
path = './fits_files/new/signal_comp_2/redu00_3evals/'
full_file_ = fits.open(path + f'toltec_commissioning_a1100_citlali.fits')
full_map = stack_fg(full_file_)[0]
full_header = full_file_[1].header
obs_names = [131947,131949,134710,134712,134714,134859,134861]
offsets_y, offsets_x, peaks_S = [], [], []
for obs in obs_names:
    fits_file_ = fits.open(path + f'toltec_commissioning_a1100_science_{str(obs)}_citlali.fits')
    map_dat = fits_file_[1].data[0][0]
    map_header = fits_file_[1].header
    print(obs)
    offset_y, offset_x, peak_S = pointing_correction(map_dat, full_map, map_header, full_header)
    offsets_y.append(offset_y)
    offsets_x.append(offset_x)
    peaks_S.append(float(f'{peak_S:.2e}'))
pointing_offsets = QTable([offsets_x, offsets_y, peaks_S], names=('x Offset', 'y Offset', 'Peak S'))
pointing_offsets.pprint_all()