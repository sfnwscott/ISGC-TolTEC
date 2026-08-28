# load packages
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import os

def point_correct(map1, header1, map2, header2, offsets1, offsets2):
    rel_offsets = offsets2 - offsets1 # in x,y
    wcs_1 = WCS(header1,naxis=2)
    wcs_2 = WCS(header2,naxis=2)
    RA, DEC = wcs_1.wcs.crval
    center_coordinate = SkyCoord(ra=RA, dec=DEC, frame='icrs', obstime="J2000", unit=(u.deg,u.deg))
    
    cy_1, cx_1 = np.round(center_coordinate.to_pixel(wcs_1))
    cy_1, cx_1 = int(cy_1), int(cx_1)
    map1 = map1[cy_1-200:cy_1+200,cx_1-200:cx_1+200]
  
    cy_2, cx_2 = np.round(center_coordinate.to_pixel(wcs_2))
    cy_2, cx_2 = int(cy_2), int(cx_2)

    map2 = map2[cy_2-160:cy_2+160,cx_2-160:cx_2+160]
    h2, w2 = map2.shape
    map1 = map1[40+rel_offsets[1]:40+h2+rel_offsets[1],40+rel_offsets[0]:40+w2+rel_offsets[0]]
    return map1, map2

def generate_psd(resid):
    fft_map = np.fft.fft2(resid)
    fft_map = np.fft.fftshift(fft_map)
    psd_2d = np.abs(fft_map)**2
    ny, nx = fft_map.shape
    fx = np.fft.fftfreq(nx, d=1) * 60
    fy = np.fft.fftfreq(ny, d=1) * 60
    fx_shifted = np.fft.fftshift(fx)
    fy_shifted = np.fft.fftshift(fy)
    FX, FY = np.meshgrid(fx_shifted, fy_shifted)
    F_radial = np.sqrt(FX**2 + FY**2)
    r = F_radial.flatten()
    P = psd_2d.flatten()
    nbins = 100
    r_bins = np.linspace(r.min(), r.max(), nbins + 1)
    r_bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    indices = np.digitize(r, r_bins)
    psd_1d = np.array([P[indices == i].mean() if np.any(indices == i) else np.nan for i in range(1, nbins + 1)])
    return r_bin_centers, psd_1d

time_start = time.time()

offset_dat = pd.read_csv('signal_table_full3eig.csv')
offsets = offset_dat.iloc[0:7][['x Offset','y Offset']].to_numpy()
psd_pairs = np.array([
    [0.25,15],
    [0.50,20],
    [1.00,50]
])
map_pairs = np.array([
    [131947,131949],
    [134710,134712],
    [134859,134861]
])

redu_neig_dict = {}
directory = os.listdir()
for path in directory:
    if not os.path.isdir(path):
        continue
    else:
        pass
    full_path = path + '/coadded/raw/toltec_commissioning_a1100_citlali.fits'
    full_fits_file = fits.open(full_path)
    full_map = full_fits_file[1].data[0][0]
    full_header = full_fits_file[1].header
    neig = full_fits_file[0].header['HIERARCH CONFIG.CLEANED.NEIG']
    redu_neig_dict[str(neig)] = path

psd_dataframe = pd.DataFrame()
for psd_pair in psd_pairs:
    ell_cutoff, neig = psd_pair
    path = redu_neig_dict[str(neig)]
    i = 0
    obs_pair_psd = []
    obs_pair_freqbins = []
    for map_pair in map_pairs:
        obs1, obs2 = map_pair
        path1 = path + '/' + str(obs1) + '/raw/toltec_commissioning_a1100_science_' + str(obs1) + '_citlali.fits'
        path2 = path + '/' + str(obs2) + '/raw/toltec_commissioning_a1100_science_' + str(obs2) + '_citlali.fits'
        fits_1 = fits.open(path1)
        fits_2 = fits.open(path2)
        map1 = fits_1[1].data[0][0]
        map2 = fits_2[1].data[0][0]
        header1 = fits_1[1].header
        header2 = fits_2[1].header
        offsets1 = offsets[0]
        offsets2 = offsets[1]
        if i == 4:
            i += 1
        map1, map2 = point_correct(map1, header1, map2, header2, offsets[i], offsets[i+1])
        i += 2
        fc_1 = map1 > 100
        fc_2 = map2 > 100
        sum_1 = np.sum(map1 * fc_1) / len((map1 * fc_1).flatten())
        sum_2 = np.sum(map2 * fc_2) / len((map2 * fc_2).flatten())
        resid = map1 - (sum_1/sum_2)*map2
        freq_bins, psd = generate_psd(resid)

        col_name = str(obs1) + '_' + str(obs2) + '_N-eig=' + str(neig)
        freqs_name = col_name + '_ell'
        p_name = col_name + 'P'
        psd_dataframe[freqs_name] = freq_bins
        psd_dataframe[p_name] = psd
psd_dataframe.to_csv('psd_data.csv')

time_end = time.time()
print(f'Time Elapsed: {time_end-time_start}s')

    
