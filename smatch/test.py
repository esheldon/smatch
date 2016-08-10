from __future__ import print_function
import sys, os
import tempfile
import unittest

import numpy


from . import smatch
from .smatch import Catalog, read_matches

def test():
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSMatch)
    unittest.TextTestRunner(verbosity=2).run(suite)

class TestSMatch(unittest.TestCase):
    def setUp(self):
        self.nside=4096

        two = 2.0/3600.
        # offset second list by fraction of 2 arcsec in dec
        # not last ones don'e match at all
        self.ra1 = numpy.array(  [200.0, 200.0, 200.0, 175.23, 21.36])
        self.dec1 = numpy.array( [24.3,          24.3,            24.3,  -28.25, -15.32])
        # make one of them big endian to check byte swapping
        self.ra2 = numpy.array(  [200.0, 200.0, 200.0, 175.23, 55.25], dtype='>f8')
        self.dec2 = numpy.array( [24.3+0.75*two, 24.3 + 0.25*two, 24.3 - 0.33*two, -28.25 + 0.58*two, 75.22])

        self.two=two

    
        self.maxmatches = [0,1,2]
        self.expected = [10,4,7]

    def testCreate(self):
        cat, ok = self.make_cat(self.two)
        self.assertTrue(ok,"creating Catalog object") 

        nside=cat.get_hpix_nside()
        self.assertEqual(nside, self.nside, "checking depth can be gotten")

        area=cat.get_hpix_area()

    def testMatch(self):

        # scalar radius
        radius=self.two
        radii = numpy.zeros(self.ra1.size) + self.two

        for rad in [radius, radii]:

            if numpy.isscalar(rad):
                rstr=' scalar radius'
            else:
                rstr=' array radius'

            cat, ok = self.make_cat(rad)
            self.assertTrue(ok,"creating Catalog object")

            if not ok:
                skipTest("cannot test result if Catalog object creation fails")

            for maxmatch, expected in zip(self.maxmatches,self.expected):
                cat.match(self.ra2, self.dec2, maxmatch=maxmatch)

                self.check_matches(cat.get_nmatches(),
                                   expected,
                                   maxmatch,
                                   rstr+' get_nmatches')
                self.check_matches(cat.get_matches().size,
                                   expected,
                                   maxmatch,
                                   rstr+' matches.size')


    def testMatch2File(self):

        # scalar radius
        radius=self.two
        radii = numpy.zeros(self.ra1.size) + self.two

        for rad in [radius, radii]:

            fname=tempfile.mktemp(prefix="testSMatch2File",suffix='.dat')

            try:
                if numpy.isscalar(rad):
                    rstr=' scalar radius'
                else:
                    rstr=' array radius'

                cat, ok = self.make_cat(rad)
                self.assertTrue(ok,"creating Catalog object")

                if not ok:
                    skipTest("cannot test result if Catalog object creation fails")

                for maxmatch, expected in zip(self.maxmatches,self.expected):
                    cat.match2file(fname, self.ra2, self.dec2, maxmatch=maxmatch)

                    self.check_matches(cat.get_nmatches(),
                                       expected,
                                       maxmatch,
                                       rstr+' get_nmatches')
                    matches=read_matches(fname)
                    self.check_matches(matches.size,
                                       expected,
                                       maxmatch,
                                       rstr+' matches.size')

            finally:
                if os.path.exists(fname):
                    os.remove(fname)

    def check_matches(self, nmatches, expected, maxmatch,extra):
        mess="expected %d matches with maxmatch=%d, got %d (%s)"
        mess = mess % (expected, maxmatch, nmatches, extra),
        self.assertEqual(
            nmatches,
            expected,
            extra,
        )

    def make_cat(self, rad):
        try:
            cat = Catalog(self.ra1, self.dec1, rad, nside=self.nside)
            ok=True
        except:
            cat = None
            ok=False

        return cat, ok

def testold(maxmatch=1, npts=100,
         depth=12, radius=0.1,
         seed=None, doplot=False):

    numpy.random.seed(seed)

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))
    cat=Catalog(ra, dec, radii)

    print(cat)
    print("initial nmatches (should be 0):",cat.nmatches)

    cat.match(ra, dec, maxmatch=maxmatch)

    print("nmatches:",cat.nmatches)

    if doplot:
        plot_points_and_radii(ra, dec, radii)

def test_file(filename, maxmatch=1, npts=100,
              depth=12, radius=0.1,
              seed=None, doplot=False):
    numpy.random.seed(seed)

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))
    cat=Catalog(ra, dec, radii)

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
                     depth=12, radius=0.1,
                     seed=None, doplot=False):
    numpy.random.seed(seed)

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))

    matches = smatch.match(ra, dec, radii, ra, dec, maxmatch=maxmatch)

    print("nmatches:",matches.size)
    if doplot:
        plot_points_and_radii(ra, dec, radii)

    return matches


def test_against_htm(maxmatch=1, npts=100,
                     depth=12, radius=0.1, seed=None, doplot=False):

    import time
    import esutil as eu

    numpy.random.seed(seed)

    depth=12

    ra = 10.0 + numpy.random.rand(npts)
    dec = 10.0 + numpy.random.rand(npts)

    #radii=numpy.zeros(ra.size) + radius
    radii=radius*(1.0 + 0.1*numpy.random.uniform(size=ra.size,low=-0.1,high=0.1))
    #radii=radius
    cat=Catalog(ra,dec,radii)
    print(cat)
    print("initial nmatches (should be 0):",cat.nmatches)

    print('Doing healpix match')
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


