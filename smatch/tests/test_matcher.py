import numpy as np
import esutil.coords

import pytest

from ..matcher import Matcher


def _gen_sphere_pts(n, seed):
    rng = np.random.RandomState(seed=seed)
    ra = rng.uniform(size=n) * 360
    dec = np.arcsin(rng.uniform(size=n, low=-1, high=1)) / np.pi * 180.0
    return ra, dec


@pytest.mark.parametrize('k', [1, 2, 3])
def test_matcher_knn(k):
    ra, dec = _gen_sphere_pts(50, 4543)
    mch = Matcher(ra, dec)
    idx, d = mch.query_knn(ra, dec, k=k, return_distances=True)

    # test via brute force
    idxs = np.arange(ra.shape[0])
    for i in range(ra.shape[0]):
        ds = esutil.coords.sphdist(ra[i], dec[i], ra, dec)
        inds = np.argsort(ds)

        if k != 1:
            assert np.allclose(d[i, :], ds[inds[:k]], atol=2e-6)
            assert np.array_equal(idx[i, :], idxs[inds[:k]])
        else:
            assert np.allclose(d[i], 0, atol=2e-6)
            assert np.array_equal(idx[i], i)


def test_matcher_knn_indices():
    ra, dec = _gen_sphere_pts(50, 4543)
    mch = Matcher(ra, dec)
    idx, i1, i2, d = mch.query_knn(ra[:-10], dec[:-10], k=1, return_indices=True)

    assert np.array_equal(i1, i2)

    # test via brute force
    for i in range(ra.shape[0]-10):
        assert np.allclose(d[i], 0, atol=2e-6)
        assert np.array_equal(idx[i], i)


@pytest.mark.parametrize('k', [1, 2, 3])
def test_matcher_knn_maxrad(k):
    ra, dec = _gen_sphere_pts(50, 4543)
    mch = Matcher(ra, dec)
    idx, d = mch.query_knn(
        ra, dec, distance_upper_bound=5e4/3600, k=k, return_distances=True
    )

    # test via brute force
    idxs = np.arange(ra.shape[0])
    for i in range(ra.shape[0]):
        ds = esutil.coords.sphdist(ra[i], dec[i], ra, dec)
        inds = np.argsort(ds)
        msk = (ds[inds] < 5e4/3600) & (np.arange(ra.shape[0]) < k)
        s = np.sum(msk)

        if k != 1:
            assert np.allclose(d[i, :s], ds[inds[msk]], atol=2e-6)
            assert np.array_equal(idx[i, :s], idxs[inds[msk]])
        else:
            assert np.allclose(d[i], 0, atol=2e-6)
            assert np.array_equal(idx[i], i)


@pytest.mark.parametrize('k', [1, 2, 3])
def test_matcher_knn_maxrad_inf(k):
    ra, dec = _gen_sphere_pts(50, 4543)
    mch = Matcher(ra, dec)
    rap, decp = _gen_sphere_pts(50, 443)
    idx, d = mch.query_knn(
        rap, decp, distance_upper_bound=1/3600, k=k, return_distances=True)
    assert not np.any(np.isfinite(d))
    assert np.all(idx == 50)


def test_matcher_radius():
    ra, dec = _gen_sphere_pts(50, 4543)
    mch = Matcher(ra, dec)

    rap, decp = _gen_sphere_pts(100, 454)
    idx = mch.query_radius(rap, decp, 6e4/3600)

    for ic in range(ra.shape[0]):
        idxc = []
        for ip in range(rap.shape[0]):
            sep = esutil.coords.sphdist(ra[ic], dec[ic], rap[ip], decp[ip])
            if sep < 6e4/3600:
                idxc.append(ip)
        assert set(idxc) == set(idx[ic])
