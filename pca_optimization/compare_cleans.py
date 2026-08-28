import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import QTable
from astropy.visualization import simple_norm
import os

path = './fits_files/new/signal_comp_2/redu00_3evals/'
obs_names = [131947,131949,134710,134712,134714]

neigs = [1,3,5,10,20,30]


fig, axs = plt.subplots(1,len(obs_names),figsize=(len(obs_names),1.5),dpi=200,constrained_layout=True)
for i,obs in enumerate(obs_names):
    fits_path = path + 'toltec_commissioning_a1100_science_' + str(obs) + '_citlali.fits'
    fits_file_ = fits.open(fits_path)
    map_dat = fits_file_[1].data[0][0]
    map_header = fits_file_[1].header
    wcs = WCS(map_header,naxis=2)
    RA, DEC = wcs.wcs.crval
    center_coordinate = SkyCoord(ra=RA, dec=DEC, frame='icrs', obstime="J2000", unit=(u.deg,u.deg))
    cy, cx = np.round(center_coordinate.to_pixel(wcs))
    cy, cx = int(cy), int(cx)
    map_dat = map_dat[cy-180:cy+180,cx-180:cx+180]
    if i == 0:
        norm = simple_norm(map_dat, min_cut=map_dat.min(), max_cut=map_dat.max())
    im = axs[i].imshow(map_dat, origin='lower', norm=norm)
    axs[i].set_xticks([])
    axs[i].set_yticks([])
    if i == len(obs_names)-1:
        cbar = fig.colorbar(im, ax=axs[i],shrink=0.65)
        cbar.ax.tick_params(labelsize=8)
    neig = fits_file_[0].header['HIERARCH CONFIG.CLEANED.NEIG']
    clean_text = r'$N_{eig}$=' + str(neig)
    prop = dict(boxstyle='round', facecolor='grey', alpha=0.5)
    axs[i].text(255,330,clean_text,color='white',fontsize=4,verticalalignment='center',bbox=prop)
plt.show()