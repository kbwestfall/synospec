.. code-block:: console

    $ wfos_etc -h
    usage: wfos_etc [-h] [--spec_file SPEC_FILE] [--spec_wave SPEC_WAVE]
                    [--spec_wave_units SPEC_WAVE_UNITS] [--spec_flux SPEC_FLUX]
                    [--spec_flux_units SPEC_FLUX_UNITS]
                    [--spec_res_indx SPEC_RES_INDX | --spec_res_value SPEC_RES_VALUE]
                    [--spec_table SPEC_TABLE] [-l EMLINE] [-m MAG]
                    [--mag_band MAG_BAND] [--mag_system MAG_SYSTEM]
                    [--sky_mag SKY_MAG] [--sky_mag_band SKY_MAG_BAND]
                    [--sky_mag_system SKY_MAG_SYSTEM] [-z REDSHIFT]
                    [-s SERSIC SERSIC SERSIC SERSIC | -u] [--refl REFL]
                    [--blue_grat BLUE_GRAT]
                    [--blue_wave BLUE_WAVE | --blue_angle BLUE_ANGLE]
                    [--blue_binning BLUE_BINNING BLUE_BINNING] [--red_grat RED_GRAT]
                    [--red_wave RED_WAVE | --red_angle RED_ANGLE]
                    [--red_binning RED_BINNING RED_BINNING]
                    [--slit SLIT SLIT SLIT SLIT SLIT] [--extract EXTRACT] [-t TIME]
                    [-f FWHM] [-a AIRMASS] [-i] [-p] [--snr_units SNR_UNITS]
                    [--sky_err SKY_ERR]
    
    WFOS Exposure-Time Calculator (v0.1)
    
    options:
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
      --spec_res_indx SPEC_RES_INDX
                            Extension or column number with the flux. (default:
                            None)
      --spec_res_value SPEC_RES_VALUE
                            Single value for the spectral resolution (R =
                            lambda/dlambda) for the full spectrum. (default: None)
      --spec_table SPEC_TABLE
                            Extension in the fits file with the binary table data.
                            (default: None)
      -l EMLINE, --emline EMLINE
                            File with emission lines to add to the spectrum.
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
      --refl REFL           Select the reflectivity curve for TMT. Must be either
                            'req' or 'goal' for the required or goal reflectivity
                            performance. (default: req)
      --blue_grat BLUE_GRAT
                            Grating to use in the blue arm. Options are: B1210,
                            B2479, B2700, B3600 (default: B1210)
      --blue_wave BLUE_WAVE
                            Central wavelength for the blue arm. If None, will use
                            the peak-efficiency wavelength. (default: None)
      --blue_angle BLUE_ANGLE
                            Grating angle for blue grating. If None, will use then
                            angle the provides the best efficiency for the on-axis
                            spectrum. (default: None)
      --blue_binning BLUE_BINNING BLUE_BINNING
                            On-chip binning for the blue grating. Order is spectral
                            then spatial. I.e., to bin 2 pixels spectrally and no
                            binning spatial, set --blue_binning 2 1 (default: [1,
                            1])
      --red_grat RED_GRAT   Grating to use in the red arm. Options are: R680, R1392,
                            R1520, R2052 (default: R680)
      --red_wave RED_WAVE   Central wavelength for the red arm. If None, will use
                            the peak-efficiency wavelength. (default: None)
      --red_angle RED_ANGLE
                            Grating angle for red grating. If None, will use then
                            angle the provides the best efficiency for the on-axis
                            spectrum. (default: None)
      --red_binning RED_BINNING RED_BINNING
                            On-chip binning for the red grating. Order is spectral
                            then spatial. I.e., to bin 2 pixels spectrally and no
                            binning spatial, set --red_binning 2 1 (default: [1, 1])
      --slit SLIT SLIT SLIT SLIT SLIT
                            Slit properties: x field center, y field center, width,
                            length, rotation. The rotation is in degrees, everything
                            else is in on-sky arcsec. The slit width is in the
                            *unrotated* frame, meaning the effective slit width for
                            a rotated slit is slit_width/cos(rotation). For the
                            field center, x is along the dispersion direction with a
                            valid range of +/- 90 arcsec, and y is in the cross-
                            dispersion direction with a valid range of +/- 249
                            arcsec. Coordinate (0,0) is on axis. (default: [0.0,
                            0.0, 0.75, 5.0, 0.0])
      --extract EXTRACT     Extraction aperture in arcsec *along the slit* centered
                            on the source. At the detector, the extraction aperture
                            is narrower by cos(slit rotation). If not provided, set
                            to the FWHM of the seeing disk; see --fwhm (default:
                            None)
      -t TIME, --time TIME  Exposure time (s) (default: 3600.0)
      -f FWHM, --fwhm FWHM  On-sky PSF FWHM (arcsec) (default: 0.65)
      -a AIRMASS, --airmass AIRMASS
                            Airmass (default: 1.0)
      -i, --ipython         After completing the setup, embed in an IPython session.
                            (default: False)
      -p, --plot            Do not provide a plot of the components of the
                            calculation. (default: True)
      --snr_units SNR_UNITS
                            The units for the S/N. Options are pixel, angstrom,
                            resolution. (default: pixel)
      --sky_err SKY_ERR     The fraction of the Poisson error in the sky incurred
                            when subtracting the sky from the observation. Set to 0
                            for a sky subtraction that adds no error to the sky-
                            subtracted spectrum; set to 1 for a sky-subtraction
                            error that is the same as the Poisson error in the sky
                            spectrum acquired during the observation. (default: 0.1)
    