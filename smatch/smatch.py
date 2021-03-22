from __future__ import print_function
from sys import stderr
import numpy as np
from . import _smatch

# area 0.013114 square degrees
NSIDE_DEFAULT=4096

def match(ra1, dec1, radius1, ra2, dec2,
          nside=NSIDE_DEFAULT, maxmatch=1,
          file=None):
    """
    match points on the sphere

    parameters
    ----------

    ra1: array
        right ascension array 1 in degrees
    dec1: array
        declination array 1 in degrees, same size as ra1
    radius: array or scalar
        search radius around each point in degrees; can be a scalar
        or same size as ra1/dec1.

    ra2: array
        right ascension array 2 in degrees
    dec2: array
        declination array 2 same size as ra2 in degrees

    nside: int, optional
        nside for the healpix layout. Default 2048

    maxmatch: int, optional
        maximum number of matches to allow per point. The closest maxmatch
        matches will be kept.  Default is 1, which implles keepin the
        closest match.  Set to <= 0 to keep all matches.

    file: string
        File in which to write matches.

    returns
    -------
    matchcat: structured array
        Structured array with fields

            i1: index in the primary Catalog
            i2: index in the second set of points
            cos(dist): cosine of the angular distance
                between the points

    If a file is sent, None is returned
    """

    cat = Catalog(ra1, dec1, radius1, nside=nside)

    cat.match(ra2, dec2, maxmatch=maxmatch, file=file)

    if file is not None:
        return None
    else:
        return cat.matches

def match_self(ra, dec, radius,
               nside=NSIDE_DEFAULT, maxmatch=1,
               file=None):
    """
    match points on the sphere.  Match the catalog to itself, 
    ignoring exact matches

    parameters
    ----------

    ra: array
        right ascension array in degrees
    dec: array
        declination array in degrees, same size as ra
    radius: array or scalar
        search radius around each point in degrees; can be a scalar
        or same size as ra1/dec1.

    nside: int, optional
        nside for the healpix layout. Default 2048

    maxmatch: int, optional
        maximum number of matches to allow per point. The closest maxmatch
        matches will be kept.  Default is 1, which implles keepin the
        closest match.  Set to <= 0 to keep all matches.

    file: string
        File in which to write matches.

    returns
    -------
    matchcat: structured array
        Structured array with fields

            i1: index in the primary Catalog
            i2: index in the second set of points
            cos(dist): cosine of the angular distance
                between the points

    If a file is sent, None is returned
    """

    cat = Catalog(ra, dec, radius, nside=nside)

    cat.match_self(maxmatch=maxmatch, file=file)

    if file is not None:
        return None
    else:
        return cat.matches


class Catalog(_smatch.Catalog):
    """
    Catalog for spacial matching on the sphere

    parameters
    ----------
    ra: array
        right ascension in degrees
    dec: array
        declination in degrees
    radius: array or scalar
        Search radius in degrees. Can be scalar or same size as ra/dec
    nside: int, optional
        nside for the healpix layout. Default 2048
    """
    def __init__(self, ra, dec, radius, nside=NSIDE_DEFAULT):

        ra,dec,radius=_get_arrays(ra,dec,radius=radius)
        self._matches = None

        #super(Catalog,self).__init__(
        #    nside, ra, dec, radius,
        #)
        super(Catalog,self).__init__(nside)
        self._ra=ra
        self._dec=dec
        self._radius=radius

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
            raise RuntimeError("no match structure found; either run "
                               "match() first or load results from a "
                               "file if you wrote matches")

        return self._matches

    def get_nmatches(self):
        """
        get the number of matches found

        This will be accurate even if matches were written to a file
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

    def match(self, ra, dec, maxmatch=1, file=None):
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
        file: filename
            Send matches to the specified file
        """

        ra,dec=_get_arrays(ra,dec)
        matching_self=0

        self._match(
            maxmatch,
            matching_self,
            ra,
            dec,
            file,
        )

    def match_self(self, maxmatch=1, file=None):
        """
        match the catalog against itself, ignoring exact
        matches

        parameters
        ----------
        maxmatch: int, optional
            maximum number of matches to allow per point. The closest maxmatch
            matches will be kept.  Default is 1, which implles keepin the
            closest match.  Set to <= 0 to keep all matches.
        """

        matching_self=1

        self._match(
            maxmatch,
            matching_self,
            self._ra,
            self._dec,
            file,
        )

    def _match(self, maxmatch, matching_self, ra, dec, file):
        """
        We keep all the logic of choosing different methods here
        """

        # make sure to store None here, since the matches are in a file
        self._matches=None

        if file is not None:
            super(Catalog, self).match2file(
                maxmatch,
                matching_self,
                self._ra,
                self._dec,
                self._radius,
                ra,
                dec,
                file,
            )

        else:
            self._matches = np.zeros(1, dtype=match_dtype)
            super(Catalog, self).match(
                maxmatch,
                matching_self,
                self._ra,
                self._dec,
                self._radius,
                ra,
                dec,
                self._matches,
            )

            nmatches = self.get_nmatches()
            assert self._matches.size == nmatches,\
                ('match count does not match: '
                 '%d in array, %d counted' % (self._matches.size, nmatches))

    def __repr__(self):
        area=self.get_hpix_area()*(180.0/np.pi)**2
        lines=[
            'smatch catalog',
            '    nside:               %d' % self.get_hpix_nside(),
            '    pixel area (sq deg): %f' % area,
            '    npoints:             %d' % self._ra.size,
        ]
        return '\n'.join(lines)

def read_matches(filename):
    """
    read matches from the indicated file

    returns
    -------
    matches: structured array
        array with fields
            i1: index in the primary Catalog
            i2: index in the second set of points
            cos(dist): cosine of the angular distance
                between the points

    """
    nmatches = _smatch._count_lines(filename)
    matches = np.zeros(nmatches, dtype=match_dtype)

    if nmatches > 0:
        _smatch._load_matches(filename, matches)
    return matches



def _get_arrays(ra, dec, radius=None):
    ra=np.array(ra, ndmin=1, dtype='f8', copy=False)
    dec=np.array(dec, ndmin=1, dtype='f8', copy=False)

    if ra.size != dec.size:
        mess="ra/dec size mismatch: %d %d"
        raise ValueError(mess % (ra.size,dec.size))

    if radius is not None:

        radarr=np.array(radius, ndmin=1, dtype='f8', copy=False)

        if radarr.size != ra.size and radarr.size != 1:
            mess=("radius has size %d but expected either "
                  "a scalar/size 1 array or array of size %d")
            raise ValueError(mess % (radarr.size,ra.size))

        return ra,dec,radarr
    else:
        return ra,dec

# data type of the match structure
match_dtype=[
    ('i1','i8'),
    ('i2','i8'),
    ('cosdist','f8'),
]
