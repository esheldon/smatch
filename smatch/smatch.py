from __future__ import print_function
from sys import stderr
import numpy
from . import _smatch


class Catalog(_smatch.Catalog):
    """
    Catalog for matching

    parameters
    ----------
    ra: array
        right ascension in degrees
    dec: array
        declination in degrees
    radius: array or scalar
        Search radius in degrees. Can be scalar or same size as ra/dec
    nside: int, optional
        nside for the healpix layout. Default 512
    """
    def __init__(self, ra, dec, radius, nside=512):

        ra,dec,radius=_get_arrays(ra,dec,radius=radius)
        self._matches = None

        super(Catalog,self).__init__(
            nside, ra, dec, radius,
        )
        self.ra=ra
        self.dec=dec
        self.radius=radius

    def get_matches(self):
        """
        get the match structure

        returns
        -------
        structure with the following structure
            i1: index in the primary Catalog
            i2: index in the second set of points
            cos(dist): cosine of the angular distance
                between the points
        """
        if self._matches is None:
            raise RuntimeError("run match() first")

        return self._matches

    def get_nmatches(self):
        """
        get the number of matches
        """
        return super(Catalog,self).get_nmatches()

    def get_hpix_nside(self):
        """
        get the nside for healpix
        """
        return super(Catalog,self).get_hpix_nside()

    def get_hpix_area(self):
        """
        get the area for healpix pixels square radians
        """
        return super(Catalog,self).get_hpix_area()


    matches=property(fget=get_matches)
    nmatches=property(fget=get_nmatches)
    hpix_nside=property(fget=get_hpix_nside)
    hpix_area=property(fget=get_hpix_nside)

    def match(self, ra, dec, maxmatch=1):
        """
        match the catalog to the second set of points

        parameters
        ----------
        ra: array
            ra to match, in degrees
        dec: array
            dec to match, in degrees
        maxmatch: int, optional
            maximum number of matches to allow per point. The closest maxmatch
            matches will be kept.  Default is 1, which implles keepin the
            closest match.  Set to <= 0 to keep all matches.
        """
        ra,dec=_get_arrays(ra,dec)
        super(Catalog, self).match(
            maxmatch,
            ra,
            dec,
        )

        nmatch=self.get_nmatches()
        matches = numpy.zeros(nmatch, dtype=_match_dtype)

        if nmatch > 0:
            super(Catalog,self)._copy_matches(matches)

        self._matches=matches

    def __repr__(self):
        area=self.get_hpix_area()*(180.0/numpy.pi)**2
        lines=[
            'smatch catalog',
            '    nside:               %d' % self.get_hpix_nside(),
            '    pixel area (sq deg): %f' % area,
            '    npoints:             %d' % self.ra.size,
        ]
        return '\n'.join(lines)

def _get_arrays(ra, dec, radius=None):
    ra=numpy.array(ra, ndmin=1, dtype='f8', copy=False)
    dec=numpy.array(dec, ndmin=1, dtype='f8', copy=False)

    if ra.size != dec.size:
        mess="ra/dec size mismatch: %d %d"
        raise ValueError(mess % (ra.size,dec.size))

    if radius is not None:

        tmprad=numpy.array(radius, ndmin=1, dtype='f8', copy=False)

        if tmprad.size == 1:
            radarr=numpy.zeros(ra.size)
            radarr[:] = tmprad[0]
        elif tmprad.size == ra.size:
            radarr=tmprad
        else:
            mess=("radius has size %d but expected either "
                  "a scalar of array of size %d")
            raise ValueError(mess % (tmprad.size,ra.size))

        return ra,dec,radarr
    else:
        return ra,dec

_match_dtype=[
    ('i1','i8'),
    ('i2','i8'),
    ('cosdist','f8'),
]
