# load packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import simple_norm
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

    map2 = map2[cy_2-150:cy_2+150,cx_2-150:cx_2+150]
    h2, w2 = map2.shape
    map1 = map1[50+rel_offsets[1]:50+h2+rel_offsets[1],50+rel_offsets[0]:50+w2+rel_offsets[0]]
    return map1, map2

time_start = time.time()

offset_dat = pd.read_csv('signal_table_fl.csv')
offsets = offset_dat.iloc[0:7][['x Offset','y Offset']].to_numpy()
psd_pairs = np.array([
    [0.10,3],
    [1,6],
    [1,10],
    [1,14],
   # [12,16],
    [12,18],
    [0.25,30],
   # [0.50,20],
   # [1.00,40],
   # [2.00,110]
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

figcol = 3
figrow = len(psd_pairs)
fig, axs = plt.subplots(figrow,figcol,figsize=(figcol*1.8,figrow*1.5),tight_layout=True,dpi=400)
path_coadd_wcs = 'redu00/coadded/raw/toltec_commissioning_a1100_citlali.fits'
fits_coadd_wcs = fits.open(path_coadd_wcs)
map_coadd = fits_coadd_wcs[1].data[0][0]
header_coadd_wcs = fits_coadd_wcs[1].header
wcs_coadd = WCS(header_coadd_wcs, naxis=2)
fig_c, axs_c = plt.subplots(figrow,1,figsize=(2,figrow*1.5),
                            tight_layout=True,dpi=400)
minmax_col = [[-25,25],[-25,25],[-25,25]]
for row,psd_pair in enumerate(psd_pairs):
    ell_cutoff, neig = psd_pair
    path = redu_neig_dict[str(int(neig))]
    i = 0
    obs_pair_psd = []
    obs_pair_freqbins = []
    for col, map_pair in enumerate(map_pairs):
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
        fc_1 = map1 > np.quantile(map1, 0.95)
        fc_2 = map2 > np.quantile(map2, 0.95)
        sum_1 = np.sum(map1 * fc_1) / len((map1 * fc_1).flatten())
        sum_2 = np.sum(map2 * fc_2) / len((map2 * fc_2).flatten())
        resid = map1 - (sum_1/sum_2)*map2
        norm = simple_norm(resid,min_cut=minmax_col[col][0],max_cut=minmax_col[col][1])
#        print(neig, map_pair, np.min(resid), np.max(resid))
        im = axs[row,col].imshow(resid,origin='lower',norm=norm)
        axs[row,col].set_xlabel(str(obs1) + ' & ' + str(obs2), fontsize=6)
        axs[row,col].set_xticks([])
        axs[row,col].set_yticks([])
        cbar = fig.colorbar(im, ax=axs[row,col],shrink=0.85)
        cbar.ax.tick_params(labelsize=5)
        clean_text = r'$N_{eig}$=' + str(neig)
        prop = dict(boxstyle='round', facecolor='grey', alpha=0.5)
        axs[row,col].text(10,15,clean_text,color='white',fontsize=4,verticalalignment='center',bbox=prop)
    path_coadd = path + '/coadded/raw/toltec_commissioning_a1100_citlali.fits'
    fits_coadd = fits.open(path_coadd)
    map_coadd = fits_coadd[1].data[0][0]
    map_coadd = map_coadd[110:440,140:470]
    norm = simple_norm(map_coadd,min_cut=-20,max_cut=80)
    im = axs_c[row].imshow(map_coadd, origin='lower', norm=norm, cmap='turbo')
    cbar = fig_c.colorbar(im, ax=axs_c[row],shrink=0.9)
    cbar.ax.tick_params(labelsize=5)
    title = r'$N_{eig}$' + f'={neig}'
    axs_c[row].set_title(title,pad=3)
    axs_c[row].set_xticks([])
    axs_c[row].set_yticks([])
   # norm = simple_norm(map_coadd, stretch='asinh', asinh_a=0.025)
   # im = axs_c[row,1].imshow(map_coadd, origin='lower', norm=norm, cmap='turbo')
   # cbar = fig_c.colorbar(im, ax=axs_c[row,1],shrink=0.85)
   # cbar.ax.tick_params(labelsize=5)
   # title = r'$N_{eig}$' + f'={neig}'
   # axs_c[row,1].set_title(title,pad=3)
   # axs_c[row,1].set_xticks([])
   # axs_c[row,1].set_yticks([])
fig.tight_layout()        
fig.savefig('noise_map_comparison_fl.pdf')
fig_c.tight_layout()        
fig_c.savefig('coadd_map_comparison_fl.pdf')
time_end = time.time()
print(f'Time Elapsed: {time_end-time_start}s')
