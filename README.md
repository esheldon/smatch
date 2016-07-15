A python code for matching points on the sphere using healpix.

examples
```python
For quick matches, use the match() function

import smatch

nside=512  # healpix nside
maxmatch=1 # return closest match
matches = smatch.match(ra1, dec2, radius, ra2, dec2,
                       nside=nside, maxmatch=maxmatch)


```
