import numpy
from astropy.io import fits

hdu = fits.open('filter_curves.fits')

# WARNING: These filter curves include the telescope and detectors!

for band in ['U', 'G', 'R', 'I', 'Z']:
    data_file = 'gunn_2001_{0}_response.db'.format(band.lower())

    numpy.savetxt(data_file, numpy.transpose([hdu[band].data['wavelength'],
                                              hdu[band].data['resnoa']]),
                  fmt=['%6.1f', '%8.4f'],
                  header='SDSS {0} filter response taken from \n'.format(band.lower()) +
                         'http://www.sdss.org/wp-content/uploads/2017/04/filter_curves.fits\n' +
                         '\nSpecifically provides the data in the \'resnoa\' column; see README' +
                         ' for more details\n' + 
                         '{0:>4s} {1:>8s}'.format('wave', 'res'))

