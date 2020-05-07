# SynOSpec Data

## Night-Sky spectral data

`sky/mkea_sky_newmoon_DEIMOS_1200_2011oct.fits.gz` is taken from
`Obs/Sky/Empirical` in [xidl](https://github.com/profxj/xidl).

units are erg/s/cm^2/angstrom/arcsec^2

----

`sky/manga/apo2maunakeasky.fits` is a rescaling of the sky data taken by
MaNGA to match this spectrum in the overlapping regions.  The spectral
resolution is lower than the DEIMOS spectrum, but the wavelength range
is larger.

mgCFrame units are 1e-17 erg/s/cm^2/angstrom/fiber

fiber area is \pi arcsec^2

----

`sky/lris_esi_skyspec_fnu.fits' is a Maunakea dark-sky spectrum provided
by Chuck Steidel.  The spectrum combines at R~1500 data from LRIS-B at
lambda <~ 4500 A with R~4000 ESI data for lambda > 4500 A.  Units are
erg/cm^2/s/Hz/arcsec^2.

----

`sky-uves.fits` and `ZenExtinct-KPNO.fits` are taken from
[desihub/specter](https://github.com/desihub/specter/blob/master/py/specter/data).


