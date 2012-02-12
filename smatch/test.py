import numpy
from ._smatch import Catalog

def test():
    ra=numpy.linspace(10,20,100)
    dec=numpy.linspace(10,20,100)

    nside=512
    maxmach=1
    #rad=numpy.array(2.0/3600., dtype='f8', ndmin=1)
    rad=0.1*numpy.random.rand(ra.size) + 2.0/3600.
    cat=Catalog(nside,maxmatch,ra,dec,rad)

    res = cat.match(ra,dec)

    print res
