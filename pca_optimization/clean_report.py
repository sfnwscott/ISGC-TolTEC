# load packages
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.modeling import models, fitting
from photutils.aperture import SkyCircularAperture, ApertureStats

def compute_map(obs_path):
    # open fits
    science_citlali = fits.open(obs_path)
    obsnum = science_citlali[0].header['OBSNUM0']

    signal, weight, sig2noise = science_citlali[1], science_citlali[2], science_citlali[6]
    signal_dat, weight_dat, sig2noise_dat = signal.data[0][0], weight.data[0][0], sig2noise.data[0][0]
    # wcs info and make aperture to clip max and min at reasonable range
    wcs_ = WCS(signal.header,naxis=2)
    RA, DEC = wcs_.wcs.crval
    center_coordinate = SkyCoord(ra=RA, dec=DEC, frame='icrs', obstime="J2000", unit=(u.deg,u.deg))
    aper = SkyCircularAperture(center_coordinate, r=0.04*u.deg)
    apstats = ApertureStats(signal_dat, aper, wcs=wcs_)
    signal_dat = np.clip(signal_dat, a_min=apstats.min-20, a_max=apstats.max+20)

    # get angular extent of plots from center
    y, x = np.indices(signal_dat.shape)
    ra, dec = wcs_.all_pix2world(x,y,0)
    coords = SkyCoord(ra*u.deg, dec*u.deg)
    dra = (coords.ra - center_coordinate.ra).to(u.arcmin)
    ddec = (coords.dec - center_coordinate.dec).to(u.arcmin)
    dra_min, dra_max = dra.value.min(), dra.value.max()
    ddec_min, ddec_max = ddec.value.min(), ddec.value.max()

    # setup subplots
    layout = [
        ["left", "top right"],
        ["left", "bottom right"],
    ]
    fig, axs = plt.subplot_mosaic(layout, figsize=(11, 5), gridspec_kw={'width_ratios': [2, 1]},
                                  constrained_layout=False, dpi=100)
    plt.subplots_adjust(hspace=0.09, wspace=0.15)
    aspect = signal_dat.shape[0] / signal_dat.shape[1]
    plot_items = [{'pos':'left','dat':signal_dat,'name':'Signal'},
                  {'pos':'top right','dat':weight_dat,'name':'Weight'},
                  {'pos':'bottom right','dat':sig2noise_dat,'name':'Signal to Noise'}]

    # plotting text with config information
    evals = science_citlali[0].header['HIERARCH CONFIG.CLEANED.NEIG']
    if science_citlali[0].header['HIERARCH CONFIG.FRUITLOOPS'] == True:
        fl_iters = science_citlali[0].header['HIERARCH FRUITLOOPS_ITER']
        fl_max = science_citlali[0].header['HIERARCH CONFIG.FRUITLOOPS.MAXITER']
    else:
        fl_iters = 0
        fl_max = science_citlali[0].header['HIERARCH CONFIG.FRUITLOOPS.MAXITER']
    mean_tau = science_citlali[0].header['MEAN_TAU']
    tau = science_citlali[0].header['HIERARCH HEADER.RADIOMETER.TAU']
    int_time = science_citlali[0].header['HIERARCH HEADER.DCS.INTEGRATIONTIME'] / 60

    texteig = 'Eigenvalues Cut: ' + str(evals)
    textother = '\n'.join((
        r'$\tau\:$:% .3f' % (tau,),
        (r'$\bar{\tau}$:% .3f' % (mean_tau,)),
        ('FruitLoops\nIters: ' + str(fl_iters) + '/' + str(fl_max))
    ))

    for item in plot_items:
        pos, dat, name = item['pos'], item['dat'], item['name']
        # show each plot type
        im = axs[pos].imshow(dat, origin='lower', cmap='inferno',
                        extent=[dra_min, dra_max, ddec_min, ddec_max],
                        aspect=aspect)
        # add text on signal plot
        if item['name'] == 'Signal':
            prop = dict(boxstyle='round', facecolor='grey', alpha=0.5)
            axs[pos].text(-4.8,4.0, texteig, color='white', fontsize=10,
            verticalalignment='top', bbox=prop)
            axs[pos].text(-4.8,-4.1, textother, color='white', fontsize=8,
            verticalalignment='bottom', bbox=prop)
        # configure axes labels
        axs[pos].set(xlabel=r"$\Delta$RA",ylabel=r"$\Delta$DEC")
        xticks = [4,3,2,1,0,-1,-2,-3,-4]
        axs[pos].set_xticks(xticks)
        axs[pos].set_xticklabels([f"{int(x)}'" for x in xticks])
        yticks = [-4,-3,-2,-1,0,1,2,3,4]
        axs[pos].set_yticks(yticks)
        axs[pos].set_yticklabels([f"{int(y)}'" for y in yticks])
        axs[pos].xaxis_inverted()
        axs[pos].set_xticklabels(axs[pos].get_xticklabels()[::-1])
        # colorbar and title
        cbar = fig.colorbar(im, ax=axs[pos], orientation='vertical', shrink=0.75, pad=0.03)
        cbar.set_label(name, rotation=270, labelpad=15, fontsize=12)
        cbar.ax.tick_params(labelsize=8)
    title = r'$\rho$-OphA 1.1mm' + '\n' + obsnum
    fig.suptitle(title, fontsize=15)

    return fig, evals


def pointing_correction(obs_path, full_path):
    obs_map = fits.open(obs_path)
    full_map = fits.open(full_path)
    obs_dat = obs_map[6].data[0][0] # 6 for S/N
    full_dat = full_map[6].data[0][0]

    # get wcs for obs and full map. clip max and min
    wcs_o = WCS(obs_map[1].header,naxis=2)
    wcs_f = WCS(full_map[1].header,naxis=2)
    RA, DEC = wcs_f.wcs.crval
    center_coordinate = SkyCoord(ra=RA, dec=DEC, frame='icrs', obstime="J2000", unit=(u.deg,u.deg))
    
    # offset to above 0 and adjust gain
    full_dat = full_dat - np.min(obs_dat)
    obs_dat = obs_dat - np.min(obs_dat)
    gain_ratio = np.max(full_dat) / np.max(obs_dat)
    obs_dat = obs_dat*gain_ratio

    # crop obs_dat -- tuned. Split into larger signal (S) map and noise (N) map
    y_min_o, x_min_o = np.int16(np.round(wcs_o.world_to_pixel(center_coordinate)) - 150)
    y_max_o, x_max_o = np.int16(np.round(wcs_o.world_to_pixel(center_coordinate)) + 150)
    obs_dat_S = obs_dat[y_min_o+30:y_max_o-80,x_min_o+75:x_max_o-35]
    obs_dat_N = obs_dat[y_min_o-50:y_max_o-20,x_min_o-10:x_max_o+20]


    fw, fh = np.shape(full_dat) # dims of full map
    # create subplots
    fig, axs = plt.subplots(2, 2, figsize=(8,7), dpi=200)
    plt.subplots_adjust(wspace=0.3,hspace=0.01)
    peak_y_f, peak_x_f = np.unravel_index(np.argmax(full_dat), full_dat.shape)
    maptypes = [obs_dat_S, obs_dat_N]
    gridlengths = [64, 16]
    stat_list = []
    # first plot maps
    for i in range(2):
        omap_dat = maptypes[i]

        # get angular extent of plot from center
        y, x = np.indices(omap_dat.shape)
        ra, dec = wcs_o.all_pix2world(x,y,0)
        coords = SkyCoord(ra*u.deg, dec*u.deg)
        dra = (coords.ra - center_coordinate.ra).to(u.arcmin)
        ddec = (coords.dec - center_coordinate.dec).to(u.arcmin)
        dra_min, dra_max = dra.value.min(), dra.value.max()
        ddec_min, ddec_max = ddec.value.min(), ddec.value.max()
        ow, oh = np.shape(omap_dat) # dims of obs map

        # location of peak signal
        peak_y_o, peak_x_o = np.unravel_index(np.argmax(omap_dat), omap_dat.shape)
    
        im = axs[i,0].imshow(omap_dat,origin='lower')
        axs[i,0].axhline(y=peak_y_o,xmin=0,xmax=350,linewidth=0.8,color='black',linestyle=':')
        axs[i,0].axvline(x=peak_x_o,ymin=0,ymax=350,linewidth=0.8,color='black',linestyle=':')

        # # configure axes labels
        axs[i,0].set(xlabel=r"$\Delta$RA",ylabel=r"$\Delta$DEC")
        labeled_xticks = np.arange(np.floor(-(ow/2)/120), np.ceil((oh/2)/120)+1)
        xticks = np.linspace(0,ow,len(labeled_xticks))
        axs[i,0].set_xticks(xticks)
        axs[i,0].set_xticklabels([f"{str(int(x))}'" for x in labeled_xticks])
        axs[i,0].set_xticklabels(axs[i,0].get_xticklabels()[::-1])
        labeled_yticks = np.arange(np.floor(-(ow/2)/120), np.ceil((oh/2)/120)+1)
        yticks = np.linspace(0,ow,len(labeled_yticks))
        axs[i,0].set_yticks(yticks)
        axs[i,0].set_yticklabels([f"{str(int(y))}'" for y in labeled_yticks])
        cbar = fig.colorbar(im, ax=axs[i,0], orientation='vertical', shrink=0.65)
        cbar.ax.tick_params(labelsize=8)
        if i == 0:
            col0_title = 'Signal Reference Map'
            col1_title = 'Cross-Correlation\nby Displacement'
        else:
            col0_title = 'Noise Reference Map'
            col1_title = 'Residual RMSE\nby Displacement'
        axs[i,0].set_title(col0_title)

        # gather rmse or cross corr
        disp_x = peak_x_f - peak_x_o
        disp_y = peak_y_f - peak_y_o
        grid_len = gridlengths[i]
        x, y = np.meshgrid(np.arange(-grid_len//2, grid_len//2),np.arange(-grid_len//2, grid_len//2))
        grid = np.stack((x,y),axis=-1).reshape((grid_len**2,2))
        statistics = []
        for disp in grid:
            obs_exp = np.zeros(full_dat.shape)
            obs_exp[disp_y+disp[0]:oh+disp_y+disp[0],disp_x+disp[1]:ow+disp_x+disp[1]] = omap_dat
            if i == 0:
                cross_corr_ = np.sum(obs_exp * full_dat) / (fw*fh) # integral f(x,y)g(x,y)dxdy / dxdy
                statistics.append(cross_corr_)
            else:
                residual = obs_exp - full_dat
                residual = residual[disp_y+disp[0]:oh+disp_y+disp[0],disp_x+disp[1]:ow+disp_x+disp[1]]
                rmse_ = np.sqrt(np.mean((residual)**2))
                statistics.append(rmse_)

        # make and plot map of each statistic
        statistics = np.array(statistics).reshape((len(statistics),1))
        stacked = np.hstack((grid, statistics))
        stat_grid = np.zeros((grid_len,grid_len))
        stat_grid[stacked[:,0].astype(int)+grid_len//2,stacked[:,1].astype(int)+grid_len//2] = stacked[:,2]
        im = axs[i,1].imshow(stat_grid,origin='lower')
        axs[i,1].set(xlabel=r"$\Delta$x",ylabel=r"$\Delta$y")
        xticks = np.linspace(0,grid_len-1,3)
        axs[i,1].set_xticks(xticks)
        axs[i,1].set_xticklabels([-grid_len//2,0,grid_len//2])
        yticks = np.linspace(0,grid_len-1,3)
        axs[i,1].set_yticks(yticks)
        axs[i,1].set_yticklabels([-grid_len//2,0,grid_len//2])
        axs[i,1].set_title(col1_title)
        if i == 0:
            name = 'Cross-Correlation'
        else:
            name = 'RMSE'
        cbar = fig.colorbar(im, ax=axs[i,1], orientation='vertical', shrink=0.65)
        cbar.ax.tick_params(labelsize=8)
        if i == 0:
            amp = np.max(stat_grid)
     
        else:
            amp = np.min(stat_grid)
        theta = 0
        center = (grid_len//2, grid_len//2)
        sigmas = (2, 1)
        
        gauss = models.Gaussian2D(amp, center[1], center[0], sigmas[1], sigmas[0], theta)
        dc = models.Const2D(np.min(stat_grid))
        fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)
        yi, xi = np.indices(stat_grid.shape)
        fit_func = fitter(gauss + dc, xi, yi, stat_grid)
        area_under_curve = 2*np.pi * fit_func.amplitude_0.value * fit_func.x_stddev_0.value * fit_func.y_stddev_0.value
        prop = dict(boxstyle='round', facecolor='grey', alpha=0.5)

        axs[i,1].text(grid_len/15,grid_len/15,str(round(area_under_curve)),color='white',
                      fontsize=8,bbox=prop)
        
    stat_list.append(area_under_curve)
    return fig


#obs_path = './fits_files/new/signal_comp_2/toltec_commissioning_a1100_science_131947_citlali.fits'
#obs_path = './fits_files/new/1eval/toltec_commissioning_a1100_science_131947_citlali.fits'
#full_path = './fits_files/new/signal_comp_2/toltec_commissioning_a1100_citlali.fits'
obs_path = sys.argv[1]
full_path = sys.argv[2]
savepath = sys.argv[3]
signal, numvals = compute_map(obs_path)
signal.savefig(savepath + f'{numvals}eval_signal.jpg')
disp = pointing_correction(obs_path, full_path)
disp.savefig(savepath + f'{numvals}eval_dispmaps.jpg')


