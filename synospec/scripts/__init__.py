
from ..util.util import all_subclasses
from . import scriptbase

from . import fobos_etc_grid
from . import fobos_etc
from . import imaging_etc
from . import wfos_etc

# Build the list of script classes
def script_classes():
    import numpy

    # Recursively collect all subclasses
    scr_c = numpy.array(list(all_subclasses(scriptbase.ScriptBase)))
    scr_n = numpy.array([c.name() for c in scr_c])
    # Construct a dictionary with the script name and class
    srt = numpy.argsort(scr_n)
    return dict([ (n,c) for n,c in zip(scr_n[srt],scr_c[srt])])

synospec_scripts = list(script_classes().keys())



