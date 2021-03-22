from __future__ import print_function
import os
import tempfile
import unittest

import numpy

from ..smatch import Catalog, read_matches, match


class TestSMatch(unittest.TestCase):
    def setUp(self):
        self.nside=4096

        two = 2.0/3600.
        # offset second list by fraction of 2 arcsec in dec
        # not last ones don'e match at all
        self.ra1  = numpy.array( [200.0, 200.0, 200.0,  175.23,  21.36])
        self.dec1 = numpy.array( [24.3,   24.3,  24.3,  -28.25, -15.32])
        # make one of them big endian to check byte swapping
        self.ra2 = numpy.array(  [        200.0,           200.0,           200.0,            175.23, 55.25], dtype='>f8')
        self.dec2 = numpy.array( [24.3+0.75*two, 24.3 + 0.25*two, 24.3 - 0.33*two, -28.25 + 0.58*two, 75.22])

        self.two=two


        self.maxmatches = [0,1,2]
        self.expected = [10,4,7]

        self.expected_self = [4, 3, 4]

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

    def testMatchScalarsWithScalarRadius(self):
        match(200, 15, self.two, self.ra1, self.dec1)

    def testMatchSelf(self):

        # scalar radius
        radius=self.two
        radii = numpy.zeros(self.ra2.size) + self.two

        for rad in [radius, radii]:

            if numpy.isscalar(rad):
                rstr=' scalar radius'
            else:
                rstr=' array radius'

            cat, ok = self.make_cat(rad, set=2)
            self.assertTrue(ok,"creating Catalog object")

            if not ok:
                skipTest("cannot test result if Catalog object creation fails")

            for maxmatch, expected in zip(self.maxmatches,self.expected_self):
                cat.match_self(maxmatch=maxmatch)

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

    def testMatchSelf2File(self):

        # scalar radius
        radius=self.two
        radii = numpy.zeros(self.ra2.size) + self.two

        for rad in [radius, radii]:
            fname=tempfile.mktemp(prefix="testSMatch2File",suffix='.dat')

            try:
                if numpy.isscalar(rad):
                    rstr=' scalar radius'
                else:
                    rstr=' array radius'

                cat, ok = self.make_cat(rad, set=2)
                self.assertTrue(ok,"creating Catalog object")

                if not ok:
                    skipTest("cannot test result if Catalog object creation fails")

                for maxmatch, expected in zip(self.maxmatches,self.expected_self):
                    cat.match_self(file=fname, maxmatch=maxmatch)

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

    def make_cat(self, rad, set=1):
        try:
            if set==1:
                ra,dec = self.ra1, self.dec1
            else:
                ra,dec = self.ra2, self.dec2
            cat = Catalog(ra, dec, rad, nside=self.nside)
            ok=True
        except:
            cat = None
            ok=False

        return cat, ok
