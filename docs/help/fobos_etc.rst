.. code-block:: console

    $ fobos_etc -h
    usage: fobos_etc [-h] [--spec_file SPEC_FILE] [--spec_wave SPEC_WAVE]
                     [--spec_wave_units SPEC_WAVE_UNITS] [--spec_flux SPEC_FLUX]
                     [--spec_flux_units SPEC_FLUX_UNITS] [--ignore_mag]
                     [--spot_fwhm SPOT_FWHM]
                     [--spec_res_indx SPEC_RES_INDX | --spec_res_value SPEC_RES_VALUE]
                     [--spec_table SPEC_TABLE] [-m MAG] [--mag_band MAG_BAND]
                     [--mag_system MAG_SYSTEM] [--sky_mag SKY_MAG]
                     [--sky_mag_band SKY_MAG_BAND] [--sky_mag_system SKY_MAG_SYSTEM]
                     [-z REDSHIFT] [-l EMLINE] [-s SERSIC SERSIC SERSIC SERSIC | -u]
                     [-t TIME] [-f FWHM] [-a AIRMASS] [-i] [-p]
                     [--snr_units SNR_UNITS] [--sky_err SKY_ERR]
    
    FOBOS Exposure-Time Calculator (v0.2)
    
    optional arguments:
      -h, --help            show this help message and exit
      --spec_file SPEC_FILE
                            A fits or ascii file with the object spectrum to use
                            (default: None)
      --spec_wave SPEC_WAVE
                            Extension or column number with the wavelengths.
                            (default: WAVE)
      --spec_wave_units SPEC_WAVE_UNITS
                            Wavelength units (default: angstrom)
      --spec_flux SPEC_FLUX
                            Extension or column number with the flux. (default:
                            FLUX)
      --spec_flux_units SPEC_FLUX_UNITS
                            Input units of the flux density. Must be interpretable
                            by astropy.units.Unit. Code assumes 1e-17 erg / (cm2 s
                            angstrom) if units are not provided. (default: None)
      --ignore_mag          Ignore the provided magnitude of the source and use the
                            provided spectrum directly. Only valid if a spectrum is
                            provided. (default: False)
      --spot_fwhm SPOT_FWHM
                            FHWM of the monochromatic spot size on the detector in
                            pixels. (default: 5.8)
      --spec_res_indx SPEC_RES_INDX
                            Extension or column number with the flux. (default:
                            None)
      --spec_res_value SPEC_RES_VALUE
                            Single value for the spectral resolution (R =
                            lambda/dlambda) for the full spectrum. (default: None)
      --spec_table SPEC_TABLE
                            Extension in the fits file with the binary table data.
                            (default: None)
      -m MAG, --mag MAG     Total apparent magnitude of the source (default: 24.0)
      --mag_band MAG_BAND   Broad-band used for the provided magnitude. Must be u,
                            g, r, i, or z. (default: g)
      --mag_system MAG_SYSTEM
                            Magnitude system. Must be either AB or Vega. (default:
                            AB)
      --sky_mag SKY_MAG     Surface brightness of the sky in mag/arcsec^2 in the
                            defined broadband. If not provided, default dark-sky
                            spectrum is used. (default: None)
      --sky_mag_band SKY_MAG_BAND
                            Broad-band used for the provided sky surface brightness.
                            Must be u, g, r, i, or z. (default: g)
      --sky_mag_system SKY_MAG_SYSTEM
                            Magnitude system. Must be either AB or Vega. (default:
                            AB)
      -z REDSHIFT, --redshift REDSHIFT
                            Redshift of the object, z (default: 0.0)
      -l EMLINE, --emline EMLINE
                            File with emission lines to add to the spectrum.
                            (default: None)
      -s SERSIC SERSIC SERSIC SERSIC, --sersic SERSIC SERSIC SERSIC SERSIC
                            Use a Sersic profile to describe the object surface-
                            brightness distribution; order must be effective radius,
                            Sersic index, ellipticity (1-b/a), position angle (deg).
                            (default: None)
      -u, --uniform         Instead of a point source or Sersic profile, assume the
                            surface brightness distribution is uniform over the
                            fiber face. If set, the provided magnitude is assumed to
                            be a surface brightness. See the MAG option. (default:
                            False)
      -t TIME, --time TIME  Exposure time (s) (default: 3600.0)
      -f FWHM, --fwhm FWHM  On-sky PSF FWHM (arcsec) (default: 0.65)
      -a AIRMASS, --airmass AIRMASS
                            Airmass (default: 1.0)
      -i, --ipython         After completing the setup, embed in an IPython session.
                            (default: False)
      -p, --plot            Provide a plot of the components of the calculation.
                            (default: False)
      --snr_units SNR_UNITS
                            The units for the S/N. Options are pixel, angstrom,
                            resolution. (default: pixel)
      --sky_err SKY_ERR     The fraction of the Poisson error in the sky incurred
                            when subtracting the sky from the observation. Set to 0
                            for a sky subtraction that adds no error to the sky-
                            subtracted spectrum; set to 1 for a sky-subtraction
                            error that is the same as the Poisson error in the sky
                            spectrum acquired during the observation. (default: 0.1)
    