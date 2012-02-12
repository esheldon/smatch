import numpy
from ._smatch import Catalog

def test(maxmatch=1):
    ra=numpy.linspace(10,10.1,100)
    dec=numpy.linspace(10,10.1,100)

    nside=512
    #rad=numpy.array(2.0/3600., dtype='f8', ndmin=1)
    rad=0.1*numpy.random.rand(ra.size) + 2.0/3600.
    cat=Catalog(nside,maxmatch,ra,dec,rad)

    mcat, minput = cat.match(ra,dec)
    #mcat, minput = cat.match(numpy.zeros(1) + 0.0,numpy.zeros(1) + 0.0)
    print 'found',mcat.size,'matches'
    #print res
