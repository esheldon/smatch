A python code for matching points on the sphere using healpix.

examples
--------

```python

#For quick matches, use the match() function

import smatch

nside=512  # healpix nside
maxmatch=1 # return closest match
matches = smatch.match(ra1, dec2, radius, ra2, dec2,
                       nside=nside, maxmatch=maxmatch)


# the match structure holds the indices of the matches
# from each data set, and the cosine of the distance
# between them
print(m.dtype.descr)
    [('i1', '<i8'), ('i2', '<i8'), ('cosdist', '<f8')]

# set maxmatch=0 to get all matches, set to a positive
# integer to return that many closest matches

all_matches = smatch.match(ra1, dec2, radius, ra2, dec2,
                           nside=nside, maxmatch=0)
some_matches = smatch.match(ra1, dec2, radius, ra2, dec2,
                            nside=nside, maxmatch=3)

# write the matches to a file, rather than keeping in memory.
# useful of the number of matches is large
smatch.match2file(ra1, dec2, radius, ra2, dec2, file="matches.dat")

# you can read them later
matches=smatch.read_matches(filename)

# To match the same data set to multiple other data sets, use
# a Catalog

cat=smatch.Catalog(ra1, dec1, radius, nside=nside)
print(cat)
smatch catalog
    nside:               512
    pixel area (sq deg): 0.013114
    npoints:             100

cat.match(ra2, dec2, maxmatch=maxmatch)
cat.match(ra3, dec3, maxmatch=maxmatch)

print("found:",cat.nmatches,"matches")
matches = cat.matches

```
