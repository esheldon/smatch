from __future__ import print_function
from sys import stderr
import numpy

from . import smatch
from .smatch import Catalog

def test(maxmatch=1, rand=True, npts=100, file=None):
    if rand:
        ra = 10.0 + 0.1*numpy.random.rand(npts)
        dec = 10.0 + 0.1*numpy.random.rand(npts)
    else:
        ra = numpy.linspace(10.0,10.1,npts)
        dec = numpy.linspace(10.0,10.1,npts)

    nside=512
    radii=0.1*numpy.random.rand(ra.size) + 2.0/3600.
    cat=Catalog(nside,ra,dec,radii)

    print('Doing match')
    cat.match(ra,dec, maxmatch=maxmatch)
    print('found',cat.get_nmatches(),'matches')
    return

    for mc, mi, d in zip(mcat, minput, dist):
        #print '%d %d %e' % (mc,mi,d*3600)
        if mc == mi:
            spc1=''
            spc2='  '
        else:
            spc1='  '
            spc2=''
        print('%s%.15e %.15e %.15e %.15e%s %2d %2d %.15e' \
                % (spc1,ra[mc],dec[mc],ra[mi],dec[mi],spc2,
                   mc,mi,d*3600))

def test_file(filename, maxmatch=1, npts=100,
              nside=512, depth=12, radius=0.1,
              seed=None, doplot=False):
    numpy.random.seed(seed)

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))
    cat=Catalog(ra, dec, radii, nside=nside)

    print(cat)
    print("initial nmatches (should be 0):",cat.nmatches)

    cat.match2file(filename, ra, dec, maxmatch=maxmatch)

    print("nmatches:",cat.nmatches)
    print("reading matches")
    matches=smatch.read_matches(filename)
    print("read:",matches.size,"matches")

    if doplot:
        plot_points_and_radii(ra, dec, radii)

    return matches

def test_convenience(maxmatch=1, npts=100,
                     nside=512, depth=12, radius=0.1,
                     seed=None, doplot=False):
    numpy.random.seed(seed)

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))

    matches = smatch.match(ra, dec, radii, ra, dec, nside=nside, maxmatch=maxmatch)

    print("nmatches:",matches.size)
    if doplot:
        plot_points_and_radii(ra, dec, radii)

    return matches


def test_against_htm(maxmatch=1, npts=100,
                     nside=512, depth=12, radius=0.1, seed=None, doplot=False):

    import time
    import esutil as eu

    numpy.random.seed(seed)

    depth=12

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    #radii=numpy.zeros(ra.size) + radius
    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))
    #radii=radius
    cat=Catalog(ra,dec,radii, nside=nside)
    print(cat)
    print("initial nmatches (should be 0):",cat.nmatches)

    print(stderr,'Doing healpix match')
    t0=time.time()
    cat.match(ra,dec, maxmatch=maxmatch)
    eu.misc.ptime(time.time()-t0)

    print('found',cat.nmatches,'matches')
    matches=cat.matches

    #eu.misc.colprint(matches['i1'],matches['i2'],matches['cosdist'],numpy.arccos(matches['cosdist'])*180.0/numpy.pi)

    print('doing htm match')
    h=eu.htm.HTM(depth)
    t0=time.time()
    m1,m2,dist=h.match(ra,dec,ra,dec,radii,maxmatch=maxmatch)
    eu.misc.ptime(time.time()-t0)

    print('found',m1.size,'matches')


    if doplot:
        plot_points_and_radii(ra, dec, radii)

def plot_points_and_radii(ra, dec, radius):
    import biggles
    plt=biggles.plot(ra,dec,visible=False)

    if numpy.isscalar(radius):
        rplot = numpy.zeros(ra.size)+radius
    else:
        rplot = radius

    circles=biggles.Circles(ra, dec, rplot, color='red')
    plt.add(circles)
    plt.show()


