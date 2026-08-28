# load packages
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import QTable
import matplotlib.pyplot as plt
import time
import os
#import signal_parse_script_fl.py

#signal_parse_script_fl

# read in offsets
signal_table = pd.read_csv('signal_table_fl.csv')

time_start = time.time()

directory = os.listdir()
obs_names = [131947,131949,134710,134712,134714,134859,134861]
neigs, obs_name, ell_cutoffs, signals_norm, signals_nonorm = [], [], [], [], []
full_path = 'redu03/coadded/raw/toltec_commissioning_a1100_citlali.fits'
full_fits_file = fits.open(full_path)
full_map = full_fits_file[1].data[0][0]
full_header = full_fits_file[1].header
offsets = signal_table.iloc[0:7][['x Offset','y Offset']].to_numpy()
for path in directory:
    if not os.path.isdir(path):
        continue
    else:
        pass
    for obs,offset in zip(obs_names,offsets):
#        path_obs = './fits_files/new/signal_comp_2/redu00_3evals/toltec_commissioning_a1100_science_' + str(obs) + '_citlali.fits'
        path_obs = path + '/' + str(obs) + '/raw/toltec_commissioning_a1100_science_' + str(obs) + '_citlali.fits'
        obs_fits_file_ = fits.open(path_obs)
        map_dat = obs_fits_file_[1].data[0][0]
        map_header = obs_fits_file_[1].header
        neig = obs_fits_file_[0].header['HIERARCH CONFIG.CLEANED.NEIG']
        
        x_offset, y_offset = offset
        wcs_o = WCS(map_header,naxis=2)
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

        full_int = full_map[100+y_offset:100+obs_h+y_offset,100+x_offset:100+obs_w+x_offset]
        row_splits_all = np.asarray(np.split(obs_map, numrows, axis=1))
        splits_all = np.asarray(np.split(row_splits_all, numcols, axis=-1))
        splits_all = np.transpose(splits_all, (2, 0, 1, 3, 4))
        row_splits_coadd = np.asarray(np.split(full_int, numrows, axis=0))
        splits_coadd = np.asarray(np.split(row_splits_coadd, numcols, axis=-1))
        CC = (i_split * splits_coadd).sum(axis=-1).sum(axis=-1)
        NCC = CC / (np.sqrt((i_split**2).sum(axis=-1).sum(axis=-1)) * np.sqrt((splits_coadd**2).sum(axis=-1).sum(axis=-1)))
        noise = 1 - NCC

    



        numrows = 10 # composite of 320
        numcols = numrows
        size = h // numrows
        row_splits_all = np.asarray(np.split(maps, numrows, axis=1))
        splits_all = np.asarray(np.split(row_splits_all, numcols, axis=-1))
        splits_all = np.transpose(splits_all, (2, 0, 1, 3, 4))
        row_splits_coadd = np.asarray(np.split(coadd, numrows, axis=0))
        splits_coadd = np.asarray(np.split(row_splits_coadd, numcols, axis=-1))
        plt.figure(figsize=(10, 6))
        clipped = np.clip(coadd, a_min=-30, a_max=30)
        plt.imshow(clipped,origin='lower')
        plt.title('Coadd')
        plt.show()
        for i in range(len(maps)):
            i_split = splits_all[i]
            CC = (i_split * splits_coadd).sum(axis=-1).sum(axis=-1)
            NCC = CC / (np.sqrt((i_split**2).sum(axis=-1).sum(axis=-1)) * np.sqrt((splits_coadd**2).sum(axis=-1).sum(axis=-1)))
            noise = 1 - NCC



        for ell_cutoff in [0.1,0.25,0.5,1,2]:
            sig_norm, sig_nonorm = map_pointing(map_dat, full_map, map_header, full_header, offset, ell_cutoff)
            signals_norm.append(sig_norm)
            signals_nonorm.append(sig_nonorm)
            neigs.append(neig)
            obs_name.append(str(obs))
            ell_cutoffs.append(ell_cutoff)

ell_table = pd.DataFrame(data=np.array([neigs,signals_norm,signals_nonorm,obs_name,ell_cutoffs]).T, columns=['N-eig','S_norm','S_nonorm','Obs','ell_cutoff'])         
ell_table.to_csv('ell_cutoff_signal_table_fl.csv')

time_end = time.time()
print(f'Time elapsed: {time_end-time_start}s')