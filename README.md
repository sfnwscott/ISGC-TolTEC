# citlali_CleanAnalysis
Scripts and notebooks used for understanding the impact of PCA cleaning on signal and noise at various spatial scales. This work was supported by the NASA Illinois Space Grant Consortium and the Northwestern Summer Undergraduate Research Grants contributing to the [TolTEC](http://toltec.astro.umass.edu/science.php) research collaboration. TolTEC is an imaging polarimeter located on the Large Millimeter Telescope (LMT) operating at 1.1, 1.4, and 2.0 mm wavelengths. 



## Atmospheric Cleaning Signal Preservation

As is the case for images captured using optical telescopes as well as radio interferometers, data from large single-dish millimeter wavelength telescopes require significant processing in order to produce science ready results with sufficiently high signal to noise ratio (SNR). This is especially the case for polarimetry, as the fraction of light from quasars, interstellar dust, and other objects that is polarized is very low. The primary source of noise for ground-based telescopes operating in the mm spectrum is the atmosphere, which absorbs orders of magnitude more millimeter light than visible light. Furthermore, water vapor emits light at millimeter wavelengths making it challenging to detect astrophysical sources over the noise. 

One method to help remedy the atmosphere's high opacity is to perform a scanning maneuver with the telescope such that the detectors view different regions of the map (image) over each 10s time chunk, essentially painting a picture over the array of highly sensitive detectors. The key assumption being leveraged is that the data associated with the largest spatial variance corresponds to the atmospheric low frequency noise, referred to as 1/f noise. The other notable source of uncertainty comes from the high frequency Gaussian, white, noise. Thus, in order to maximize SNR, we must remove the large uniform structures from the map.

This is where the widespread statistical technique of Principal Component Analysis (PCA) comes in. PCA is most often used to reduce the dimensionality of a dataset by identifying the components with the most spread and subtracting away the directions with the high or low variance. The atmospheric cleaning stage of the data reduction pipeline known as `citlali` takes time ordered data, decomposes it into its principal components (its eigenvectors) and removes the most N components with the highest spatial variance, with the hope that little of the true signal is stripped away in the process. While this works great for point sources like the quasar 3C286, it may present challenges when the source itself is extended such as the star forming region Rho-Oph A spanning ~5 arcminutes. See [Rho Ophiuchi A](#rho-ophiuci-a) for more information. My goal was to better understand the impact that the atmospheric cleaning algorithm has on both the signal and noise of a map, so that we can optimize the configuration of the reduction pipeline to preserve as much astrophysical signal as possible. While in theory we could easily remove nearly all atmospheric 1/f noise, the result would be a deconstructed mess of an image lacking extended spatial structure, the primary benefit contribution of a single dish instrument like the LMT. 

### My Approach

I ran hundreds of reductions on maps varying a number of parameters, most importantly being the number of eigenvalues, or components, removed during cleaning. Determining signal presented a difficulty because there was no "truth" image to compare to, so I used a lightly cleaned map composed of multiple others to provide a high signal reference map to be compared with others. Then, I implemented a manual cross-correlation algorithm to precisely align images and get a baseline for how well they agree. As will be explained later in [#], maps varied in alignment and could not be directly superposed. Signal could be determined using cross correlation and normalized for comparison.

Because PCA cleaning removes the largest structures first, it's important to isolate these scales to see how they vary individually. I used a 2D Fourier Transform to analyze and manipulate the images. For each spatial frequency, 0.25, 0.5, 1.0, and 2.0 $\text{arcmin}^{-1}$, I removed all scales smaller than that from the maps before calculating the signal. That way, we could examine the decay in signal corresponding to certain sized shapes and larger. 

For noise, maps would similarly be precisely aligned before being subtracted from one another. Theoretically, observations taken at different times should have the same signal (to a scaling factor) with separate noise elements. Signal-less maps should include only noise, which can be characterized according to its spatial frequency by producing a Power Spectral Density plot, a tool for examining how much of each spatial scale is present. The right side corresponds to the highest frequencies like Gaussian noise while the left shows the low frequency large scale noise. Unfortunately, the isolation of true spatial noise is exceptionally difficult and many of the methods I used to remove the source such as isophotes, PCA, and z-score masking did not produce trustworthy enough results. While insights were certainly gained from examining the relationships present in the sky noise, I refocused my attention more towards determining how signal is attenuated.




## Pointing

Observations are broken into separate maps taken in 20 minute exposures that are coadded together to boost SNR. In order to do this, precise positioning information must be found to ensure the structure remains intact and the averaging process doesn't blur it. The primary method for doing this is calibration observations on sources with well-known positions to tell the 50 m telescope exactly where it is before viewing science targets. The pipeline also takes various steps to line up maps, but it isn't perfect.

Using a similar cross-correlation methodology as before to isolate precise pointing vector offsets, I improved the pointing correction operation of citlali, recovering 11% more flux than prior coadded maps. 

## Feathering

Every astronomical instrument fulfills a different need. Larger diameter telescopes are capable of imaging at higher angular resolution, allowing the study of distant protoplanetary systems, black holes, ancient receding galaxies, and more. In the millimeter and radio wavelengths, the method of interferometry has been developed, utilizing dozens of individual telescopes to synthesize an image of resolutions as fine as 5 milliarcseconds. Astronomers are able to probe distant radio sources and learn about star forming regions invisible to optical telescopes. A cost of this is spatial scale, which is where single dish observatories like the LMT excel. The natural question emerges of whether these strengths can be combined; they can, using a technique called feathering.

Feathering involves transforming single-dish images and interferometer data by first interpolating the two onto a common grid. They are then transformed to the spatial frequency domain where each is weighted according to its resolution to ensure proper preservation of spatial scale. Finally, they're combined and converted back to the spatial plane. It's a powerful tool, but can make calibrating scales challenging, so my implementation was largely a proof-of-concept rather than a practical result. It gave me the opportunity to learn about interferometry and use the CASA software, a powerful toolkit developed for the Atacama Large Millimeter Array. Unfortunately, results cannot be shown here due to 

## Data Availability

The time ordered data and reduction outputs referenced in this analysis are not included in this repository, as this work is part of an ongoing study that has not yet been published. In the meantime, this repository focuses on the analysis scripts and methodology used to produce the results described above.

## Rho Ophiuchi A

## Rho Ophiuchi A

Rho Ophiuchi A is a nearby star-forming region located within the larger Rho Ophiuchi molecular cloud complex, roughly 140 parsecs from Earth. Its proximity and youth make it one of the most heavily studied sites of active, low-mass star formation, hosting a dense cluster of young stellar objects still embedded in their natal gas and dust.

Dust grains within the cloud tend to align with the local magnetic field, causing the thermal emission they radiate to be polarized. By mapping this polarized emission at millimeter wavelengths, instruments like TolTEC can reconstruct the plane-of-sky structure of the magnetic field threading the cloud. This is important because magnetic fields are thought to play a significant role in regulating how dense cores collapse to form stars, competing with turbulence and gravity to set the efficiency and rate of star formation. Because Rho Ophiuchi A is extended (spanning roughly 5 arcminutes), it also serves as a useful test case for evaluating how atmospheric cleaning affects extended polarized structure rather than just compact point sources.