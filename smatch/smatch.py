from __future__ import print_function
from sys import stderr
import numpy
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
        nside for the healpix layout. Default 512

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
            i2: index in the second set of points (sent as arguments
                to match2file)
            cos(dist): cosine of the angular distance
                between the points

    If a file is sent, None is returned
    """

    # add early check  so we don't fail after building the
    # entire catalog, which could take a while

    ra2,dec2 = _get_arrays(ra2,dec2)


    cat = Catalog(ra1, dec1, radius1, nside=nside)

    cat.match(ra2, dec2, maxmatch=maxmatch, file=file)

    if file is not None:
        return None
    else:
        return cat.get_matches()

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
        nside for the healpix layout. Default 512

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
            i2: index in the second set of points (sent as arguments
                to match2file)
            cos(dist): cosine of the angular distance
                between the points

    If a file is sent, None is returned
    """

    cat = Catalog(ra, dec, radius, nside=nside)

    cat.match_self(maxmatch=maxmatch, file=file)

    if file is not None:
        return None
    else:
        return cat.get_matches()


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
        nside for the healpix layout. Default 512
    """
    def __init__(self, ra, dec, radius, nside=NSIDE_DEFAULT):

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
            i2: index in the second set of points (sent to match or
                match2file)
            cos(dist): cosine of the angular distance
                between the points
        """
        if self._matches is None:
            raise RuntimeError("no match structure found; either run "
                               "match() first or load results from a "
                               "file if you ran match2file()")

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
            self.ra,
            self.dec,
            file,
        )

    def match2file(self, filename, ra, dec, maxmatch=1):
        """
        deprecated, use match(..., file=filename)

        match the catalog to the second set of points, writing
        results to a file


        parameters
        ----------
        filename: string
            File in which to write the matches
        ra: array
            ra to match, in degrees
        dec: array
            dec to match, in degrees
        maxmatch: int, optional
            maximum number of matches to allow per point. The closest maxmatch
            matches will be kept.  Default is 1, which implles keepin the
            closest match.  Set to <= 0 to keep all matches.
        """

        self.match(ra, dec, maxmatch=maxmatch, file=filename)

    def _match(self, maxmatch, matching_self, ra, dec, file):
        """
        We keep all the logic of choosing different methods here
        """

        if file is not None:
            super(Catalog, self).match2file(
                maxmatch,
                matching_self,
                ra,
                dec,
                file,
            )

            # make sure to store None here, since the matches are in a file
            self._matches=None

        else:
            super(Catalog, self).match(
                maxmatch,
                matching_self,
                ra,
                dec,
            )

            nmatches = self.get_nmatches()
            matches  = numpy.zeros(nmatches, dtype=match_dtype)

            if nmatches > 0:
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

def read_matches(filename):
    """
    read matches from the indicated file

    returns
    -------
    matches: structured array
        array with fields
            i1: index in the primary Catalog
            i2: index in the second set of points (sent as arguments
                to match2file)
            cos(dist): cosine of the angular distance
                between the points

    """
    nmatches = _smatch._count_lines(filename)
    matches = numpy.zeros(nmatches, dtype=match_dtype)

    if nmatches > 0:
        _smatch._load_matches(filename, matches)
    return matches



def _get_arrays(ra, dec, radius=None):
    ra=numpy.array(ra, ndmin=1, dtype='f8', copy=False)
    dec=numpy.array(dec, ndmin=1, dtype='f8', copy=False)

    if ra.size != dec.size:
        mess="ra/dec size mismatch: %d %d"
        raise ValueError(mess % (ra.size,dec.size))

    if radius is not None:

        tmprad=numpy.array(radius, ndmin=1, dtype='f8', copy=False)

        if tmprad.size == ra.size:
            radarr=tmprad
        elif tmprad.size == 1:
            radarr=numpy.zeros(ra.size)
            radarr[:] = tmprad[0]
        else:
            mess=("radius has size %d but expected either "
                  "a scalar of array of size %d")
            raise ValueError(mess % (tmprad.size,ra.size))

        return ra,dec,radarr
    else:
        return ra,dec

# data type of the match structure
match_dtype=[
    ('i1','i8'),
    ('i2','i8'),
    ('cosdist','f8'),
]
