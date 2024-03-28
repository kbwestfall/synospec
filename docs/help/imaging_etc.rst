.. code-block:: console

    $ imaging_etc -h
    usage: imaging_etc [-h] [-m MAG] [--mag_band MAG_BAND] [--mag_system MAG_SYSTEM]
                       [-s SERSIC SERSIC SERSIC SERSIC] [--binning BINNING]
                       [-t TIME] [-n NEXP] [-f FWHM] [--beta BETA] [-a AIRMASS]
    
    Imaging Exposure-Time Calculator
    
    options:
      -h, --help            show this help message and exit
      -m MAG, --mag MAG     Total apparent magnitude of the source (default: 24.0)
      --mag_band MAG_BAND   Broad-band used for the provided magnitude. Must be u,
                            g, r, i, or z. (default: g)
      --mag_system MAG_SYSTEM
                            Magnitude system. Must be either AB or Vega. (default:
                            AB)
      -s SERSIC SERSIC SERSIC SERSIC, --sersic SERSIC SERSIC SERSIC SERSIC
                            Use a Sersic profile to describe the object surface-
                            brightness distribution; order must be effective radius,
                            Sersic index, ellipticity (1-b/a), position angle (deg).
                            If None, assume a point source. (default: None)
      --binning BINNING     Number of binned pixels along one dimension; e.g.,
                            setting "-binning 4" means 4x4 binning. (default: 8)
      -t TIME, --time TIME  Total exposure time for all combined exposures (s)
                            (default: 30.0)
      -n NEXP, --nexp NEXP  Number of exposures (default: 2)
      -f FWHM, --fwhm FWHM  On-sky PSF FWHM (arcsec) (default: 0.65)
      --beta BETA           Beta for Moffat seeing function. If None, use a
                            Gaussian. (default: None)
      -a AIRMASS, --airmass AIRMASS
                            Airmass (default: 1.0)
    