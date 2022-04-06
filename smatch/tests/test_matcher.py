import numpy as np

import pytest

from ..matcher import Matcher, sphdist


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
        ds = sphdist(ra[i], dec[i], ra, dec)
        inds = np.argsort(ds)

        if k != 1:
            assert np.allclose(d[i, :], ds[inds[:k]], atol=3e-6)
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

    idx, i1, i2, d = mch.query_knn([10], [5], k=1, return_indices=True)
    ds = sphdist([10], [5], ra, dec)
    tind = np.argmin(ds)
    assert idx == [[tind]]
    assert i1[0] == tind
    assert i2[0] == 0
    assert i1.shape == (1,)
    assert i2.shape == (1,)
    assert np.allclose(ds[tind], d[0])


def test_matcher_knn_indices_nomatch():
    ra = np.arange(10, dtype=np.float64)
    dec = np.arange(10, dtype=np.float64)

    mch = Matcher([30.0], [30.0])

    idx, i1, i2, d = mch.query_knn(ra, dec, k=1, distance_upper_bound=1.0, return_indices=True)

    assert len(i1) == 0
    assert len(i2) == 0
    assert np.all(~np.isfinite(d))


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
        ds = sphdist(ra[i], dec[i], ra, dec)
        inds = np.argsort(ds)
        msk = (ds[inds] < 5e4/3600) & (np.arange(ra.shape[0]) < k)
        s = np.sum(msk)

        if k != 1:
            assert np.allclose(d[i, :s], ds[inds[msk]], atol=3e-6)
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
            sep = sphdist(ra[ic], dec[ic], rap[ip], decp[ip])
            if sep < 6e4/3600:
                idxc.append(ip)
        assert set(idxc) == set(idx[ic])


def test_matcher_radius_indices():
    ra, dec = _gen_sphere_pts(50, 4543)
    mch = Matcher(ra, dec)

    rap = ra[::-1] + 0.1
    decp = dec[::-1] + 0.1

    idx, i1, i2, d = mch.query_radius(rap, decp, 0.2, return_indices=True)

    assert np.max(d) < 0.2
    dist_match = sphdist(ra[i1], dec[i1], rap[i2], decp[i2])
    assert np.max(dist_match) < 0.2


def test_match_radius_nomatch():
    ra = np.arange(10, dtype=np.float64)
    dec = np.arange(10, dtype=np.float64)

    mch = Matcher([30.0], [30.0])

    idx = mch.query_radius(ra, dec, 0.2)

    assert len(idx) == 1
    assert len(idx[0]) == 0


def test_match_radius_indices_nomatch():
    ra = np.arange(10, dtype=np.float64)
    dec = np.arange(10, dtype=np.float64)

    mch = Matcher([30.0], [30.0])

    idx, i1, i2, d = mch.query_radius(ra, dec, 0.2, return_indices=True)

    assert len(idx) == 1
    assert len(idx[0]) == 0
    assert len(i1) == 0
    assert len(i2) == 0
    assert len(d) == 0


def test_sphdist():
    d = sphdist(10, 20, 10, 21)
    assert np.shape(d) == tuple()
    assert np.allclose(d, 1, atol=1e-6)

    d = sphdist(10, 20, [10, 10], [21, 21.5])
    assert np.shape(d) == (2,)
    assert np.allclose(d, [1, 1.5], atol=1e-6)

    d = sphdist(10, 0, [10, 11], [0, 0])
    assert np.shape(d) == (2,)
    assert np.allclose(d, [0, 1], atol=1e-6)

    d = sphdist([10, 10], [0, 0], [10, 11], [0, 0])
    assert np.shape(d) == (2,)
    assert np.allclose(d, [0, 1], atol=1e-6)

    with pytest.raises(ValueError):
        sphdist([10], [0, 0], [10, 11], [0, 0])

    with pytest.raises(ValueError):
        sphdist(10, [0, 0], [10, 11], [0, 0])

    with pytest.raises(ValueError):
        sphdist([10, 10], [0, 0], [10, 11], 0)

    with pytest.raises(ValueError):
        sphdist([10, 10], [0, 0], [10, 11], [0])

    d = sphdist(10, 20, [[10, 10]], [[21, 21.5]])
    assert np.shape(d) == (1, 2)
    assert np.allclose(d, [[1, 1.5]], atol=1e-6)

    d = sphdist([10], [20], [[10, 10]], [[21, 21.5]])
    assert np.shape(d) == (1, 2)
    assert np.allclose(d, [[1, 1.5]], atol=1e-6)

    d = sphdist([[10, 10]], [[20, 20]], [[10, 10]], [[21, 21.5]])
    assert np.shape(d) == (1, 2)
    assert np.allclose(d, [[1, 1.5]], atol=1e-6)

    d = sphdist([[10, 10], [10, 10]], [[20, 20], [21, 21]], [[10, 10]], [[21, 21.5]])
    assert np.shape(d) == (2, 2)
    assert np.allclose(d, [[1, 1.5], [0, 0.5]], atol=2e-6)
