.. code-block:: console

    $ fobos_etc_grid -h
    usage: fobos_etc_grid [-h] [-w WAVELENGTHS WAVELENGTHS WAVELENGTHS]
                          [-m MAG MAG MAG] [-t TIME TIME TIME] [-f FWHM]
                          [-a AIRMASS] [--snr_units SNR_UNITS]
    
    Grid sampling of FOBOS ETC
    
    optional arguments:
      -h, --help            show this help message and exit
      -w WAVELENGTHS WAVELENGTHS WAVELENGTHS, --wavelengths WAVELENGTHS WAVELENGTHS WAVELENGTHS
                            Wavelength grid: start wave, approx end wave,
                            logarithmic step. (default: [3100, 10000, 4e-05])
      -m MAG MAG MAG, --mag MAG MAG MAG
                            Object apparent AB magnitude; assumed spectrum has a
                            constant AB otherwise, it is assumed to be the total
                            magnitude of the source. (default: [21.0, 25.0, 0.5])
      -t TIME TIME TIME, --time TIME TIME TIME
                            Exposure time (s) (default: [300.0, 3600.0, 300.0])
      -f FWHM, --fwhm FWHM  On-sky PSF FWHM (arcsec) (default: 0.65)
      -a AIRMASS, --airmass AIRMASS
                            Airmass (default: 1.0)
      --snr_units SNR_UNITS
                            The units for the S/N. Options are pixel, angstrom,
                            resolution. (default: pixel)
    