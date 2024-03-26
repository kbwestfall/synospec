
from pathlib import Path

from IPython import embed

import numpy
from matplotlib import pyplot

def main():
    
    files = ['super_u_asbuilt.txt', 'super_g_asbuilt.txt', 'super_r_asbuilt.txt',
             'super_i_asbuilt.txt', 'super_z_asbuilt.txt']

    for f in files:
        ifile = Path(f).resolve()
        ofile = ifile.parent / f'{ifile.with_suffix("")}_convert.txt'

        db = numpy.genfromtxt(ifile)
        srt = numpy.argsort(db[:,0])
        db[:,0] = db[srt,0]
        db[:,1] = db[srt,1]/100
        numpy.savetxt(ofile, db, fmt=['%7.1f', '%18.15f'], header=f'{"WAVE":>5} {"EFF":>18}')
        pyplot.plot(db[:,0], db[:,1])

    pyplot.show()

if __name__ == '__main__':
    main()

