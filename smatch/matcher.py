import functools
import operator

from scipy.spatial import cKDTree
import numpy as np


def _lonlat2vec(lon, lat):
    lonr = np.deg2rad(lon)
    latr = np.deg2rad(lat)
    coslat = np.cos(latr)
    return np.stack([
        np.atleast_1d(np.cos(lonr) * coslat),
        np.atleast_1d(np.sin(lonr) * coslat),
        np.atleast_1d(np.sin(latr)),
    ], axis=-1)


def sphdist(lon1, lat1, lon2, lat2):
    """Return the great circle arc distance between two sets of points.

    Note that (lon1, lat1) and (lon2, lat2) must be broadcastable.

    Parameters
    ----------
    lon1 : float or array-like
        The right ascension of the first set in degrees.
    lat1 : float or array-like
        The declination of the first set in degrees.
    lon2 : float or array-like
        The right ascension of the second set in degrees.
    lat2 : float or array-like
        The declination of the second set in degrees.

    Returns
    -------
    d : float or array-like
        The great circle arc distance between the sets in degrees.
    """
    if np.shape(lon1) != np.shape(lat1):
        raise ValueError(
            "lon1 and lat1 must be the same shape for sphdist: lon1=%s lat1=%s" % (
                np.shape(lon1), np.shape(lat1)
            )
        )

    if np.shape(lon2) != np.shape(lat2):
        raise ValueError(
            "lon2 and lat2 must be the same shape for sphdist: lon2=%s lat2=%s" % (
                np.shape(lon1), np.shape(lat1)
            )
        )

    if (
        np.shape(lon1) == tuple()
        and np.shape(lat1) == tuple()
        and np.shape(lon2) == tuple()
        and np.shape(lat2) == tuple()
    ):
        is_scalar = True
    else:
        is_scalar = False

    vec1 = _lonlat2vec(lon1, lat1)
    vec2 = _lonlat2vec(lon2, lat2)

    cosd = (
        vec1[..., 0] * vec2[..., 0]
        + vec1[..., 1] * vec2[..., 1]
        + vec1[..., 2] * vec2[..., 2]
    )

    np.clip(cosd, -1, 1, out=cosd)

    d = np.rad2deg(np.arccos(cosd))

    if is_scalar:
        return np.ravel(d)[0]
    else:
        return d


class Matcher(object):
    """A class to match catalogs with cKDTree.

    A collaboration of Alex Drlica-Wagner, Matthew Becker, & Eli Rykoff.

    Parameters
    ----------
    lon : array-like
        The longitude in degrees.
    lat : array-like
        The latitude in degrees.
    """
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat
        coords = _lonlat2vec(lon, lat)
        # The tree in the match does not need to be balanced, and
        # turning this off yields significantly faster runtime.
        self.tree = cKDTree(coords, compact_nodes=False, balanced_tree=False)

    def query_knn(
        self, lon, lat, k=1, distance_upper_bound=None, eps=0,
        return_indices=False, return_distances=False,
    ):
        """Find the `k` nearest neighbors of each point in (lon, lat) in
        the points held by the matcher.

        Parameters
        ----------
        lon : array-like, float
            The longitude in degrees.
        lat : array-like, float
            The latitude in degrees.
        k : `int`, optional
            The number of nearest neighbors to find.
        distance_upper_bound : `float`, optional
            The maximum allowed distance in degrees for a nearest neighbor.
            Default of None results in no upper bound on the distance.
        eps : `float`, optional
            If non-zero, the set of returned points are correct to within a
            fraction precision of `eps` being the actual knn
        return_indices : `bool`, optional
            Return tuple of (idx, i1, i2, d) instead of idx.
            Only supported with k=1.
        return_distances : `bool`, optional
            If True, the array of distances is returned as the last return value.
            Implied and thus ignored if `return_indices` is True.

        Returns
        -------
        idx : array-like, int
            Array of indices.  Same shape as input array with axis of
            dimension `k` added to the end.  If `k=1` then this last
            dimension is squeezed out.
        i1 : array-like, int
            Array of indices for matcher lon/lat.
            Returned if return_indices is True.
        i2 : array-like, int
            Array of indices for query lon/lat.
            Returned if return_indices is True.
        d : array-like, float
            Array of distances in degrees. Same shape as input array with axis
            of dimension `k` added to the end. If `k=1`, then this last dimension
            is squeezed out.
        """
        if distance_upper_bound is not None:
            maxd = 2*np.sin(np.deg2rad(distance_upper_bound)/2.)
        else:
            maxd = np.inf

        if k != 1 and return_indices:
            raise NotImplementedError("Indices are only returned for 1-1 matches")

        coords = _lonlat2vec(lon, lat)
        d, idx = self.tree.query(coords, k=k, p=2, distance_upper_bound=maxd, eps=eps)

        if return_distances or return_indices:
            d /= 2
            np.arcsin(d, out=d, where=np.isfinite(d))
            d = np.rad2deg(2*d)

        if return_indices:
            i2, = np.where(np.isfinite(d))
            i1 = idx[i2]
            return idx, i1, i2, d
        else:
            if return_distances:
                return idx, d
            else:
                return idx

    def query_radius(self, lon, lat, radius, eps=0.0, return_indices=False):
        """Find all points in (lon, lat) that are within `radius` of the
        points in the matcher.

        Parameters
        ----------
        lon : array-like
            The longitude in degrees.
        lat : array-like
            The latitude in degrees.
        radius : `float`
            The match radius in degrees.
        eps : `float`, optional
            If non-zero, the set of returned points are correct to within a
            fraction precision of `eps` being closer than `radius`.
        return_indices : `bool`, optional
            Return tuple of (idx, i1, i2, distance) instead of just idx.

        Returns
        -------
        idx : `list` [`list` [`int`]]
            Each row in idx corresponds to each position in matcher lon/lat.
            The indices in the row correspond to the indices in query lon/lat.
        i1 : array-like, int
            Array of indices for matcher lon/lat.
            Returned if return_indices is True.
        i2 : array-like, int
            Array of indices for query lon/lat.
            Returned if return_indices is True.
        distance : array-like, float
            Array of distance (degrees) for each match pair.
            Returned if return_indices is True.
        """
        coords = _lonlat2vec(lon, lat)
        # The second tree in the match does not need to be balanced, and
        # turning this off yields significantly faster runtime.
        qtree = cKDTree(coords, compact_nodes=False, balanced_tree=False)
        angle = 2.0*np.sin(np.deg2rad(radius)/2.0)
        idx = self.tree.query_ball_tree(qtree, angle, eps=eps)

        if return_indices:
            n_match_per_obj = np.array([len(row) for row in idx])
            i1 = np.repeat(np.arange(len(idx)), n_match_per_obj)
            i2 = np.array(functools.reduce(operator.iconcat, idx, []))
            if len(i1) > 0:
                ds = sphdist(
                    self.lon[i1], self.lat[i1],
                    lon[i2], lat[i2],
                )
            else:
                ds = np.zeros(0)
            return idx, i1, i2, ds
        else:
            return idx

    def query_self(self, radius, min_match=1, eps=0.0, return_indices=False):
        """Match the list of lon/lat to itself.

        Parameters
        ----------
        radius : `float`
            The match radius in degrees.
        min_match : `int`, optional
            Minimum number of matches to count as a match.
            If min_match=1 then all positions will be returned since every
            position will match at least to itself.
        eps : `float`, optional
            If non-zero, the set of returned points are correct to within a
            fraction precision of `eps` being closer than `radius`.
        return_indices : `bool`, optional
            Return tuple of (idx, i1, i2, distance) instead of just idx.

        Returns
        -------
        idx : `list` [`list` [`int`]]
            Each row in idx corresponds to each position in matcher lon/lat.
            The indices in the row correspond to the indices in query lon/lat.
        i1 : array-like
            Array of indices for matcher lon/lat.
            Returned if return_indices is True.
        i2 : array-like
            Array of indices for query lon/lat.
            Returned if return_indices is True.
        distance : array-like
            Array of distance (degrees) for each match pair.
            Returned if return_indices is True.
        """
        angle = 2.0*np.sin(np.deg2rad(radius)/2.0)
        idx = self.tree.query_ball_tree(self.tree, angle, eps=eps)
        if return_indices:
            n_match_per_obj = np.array([len(row) for row in idx])
            i1 = np.repeat(np.arange(len(idx)), n_match_per_obj)
            i2 = np.array(functools.reduce(operator.iconcat, idx, []))
            if len(i1) > 0:
                ds = sphdist(
                    self.lon[i1], self.lat[i1],
                    self.lon[i2], self.lat[i2],
                )
            else:
                ds = np.zeros(0)
            return idx, i1, i2, ds
        else:
            return idx

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        # Clear the memory from the tree
        del self.tree
