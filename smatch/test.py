from sys import stderr
import numpy
from ._smatch import Catalog

def test_norad(maxmatch=1, rand=True, npts=100, verbose=False):
    dorad=0
    if rand:
        ra = 10.0 + 0.1*numpy.random.rand(npts)
        dec = 10.0 + 0.1*numpy.random.rand(npts)
    else:
        ra = numpy.linspace(10.0,10.1,npts)
        dec = numpy.linspace(10.0,10.1,npts)

    nside=512
    rad=0.1*numpy.random.rand(ra.size) + 2.0/3600.
    cat=Catalog(nside,maxmatch,ra,dec,rad)

    print >>stderr,'Doing match'
    mcat, minput = cat.match(ra,dec,dorad)
    print 'found',mcat.size,'matches'

    if not verbose:
        del cat
        del mcat
        del minput
        return
    for mc, mi in zip(mcat, minput):
        print '%2d %2d' % (mc,mi)

    del cat
    del mcat
    del minput

def test_rad(maxmatch=1, rand=True, npts=100, verbose=False):
    dorad=1
    if rand:
        ra = 10.0 + 0.1*numpy.random.rand(npts)
        dec = 10.0 + 0.1*numpy.random.rand(npts)
    else:
        ra = numpy.linspace(10.0,10.1,npts)
        dec = numpy.linspace(10.0,10.1,npts)

    nside=512
    rad=0.1*numpy.random.rand(ra.size) + 2.0/3600.
    cat=Catalog(nside,maxmatch,ra,dec,rad)

    print >>stderr,'Doing match'
    mcat, minput, dist = cat.match(ra,dec,dorad)
    print 'found',mcat.size,'matches'

    if not verbose:
        del cat
        del mcat
        del minput
        del dist
        return
    for mc, mi, d in zip(mcat, minput, dist):
        #print '%d %d %e' % (mc,mi,d*3600)
        if mc == mi:
            spc1=''
            spc2='  '
        else:
            spc1='  '
            spc2=''
        print '%s%.15e %.15e %.15e %.15e%s %2d %2d %.15e' \
                % (spc1,ra[mc],dec[mc],ra[mi],dec[mi],spc2,
                   mc,mi,d*3600)

    del cat
    del mcat
    del minput
    del dist


def test_against_htm(maxmatch=1, rand=True, npts=100, verbose=False, nside=512, depth=12):
    import time
    import esutil as eu

    nside=1024
    depth=12

    dorad=0
    if rand:
        ra = 10.0 + 0.1*numpy.random.rand(npts)
        dec = 10.0 + 0.1*numpy.random.rand(npts)
    else:
        ra = numpy.linspace(10.0,10.1,npts)
        dec = numpy.linspace(10.0,10.1,npts)

    rad=0.1*numpy.random.rand(ra.size) + 2.0/3600.
    cat=Catalog(nside,maxmatch,ra,dec,rad)

    print >>stderr,'Doing healpix match'
    t0=time.time()
    mcat, minput = cat.match(ra,dec,dorad)
    eu.misc.ptime(time.time()-t0)

    print 'found',mcat.size,'matches'

    del cat
    del mcat
    del minput


    print >>stderr,'doing htm match'
    h=eu.htm.HTM(depth)
    t0=time.time()
    m1,m2,dist=h.match(ra,dec,ra,dec,rad,maxmatch=maxmatch)
    eu.misc.ptime(time.time()-t0)

    print 'found',m1.size,'matches'
