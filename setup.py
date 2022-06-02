from setuptools import setup, Extension, find_packages
import numpy


ext = Extension(
    "smatch._smatch",
    ["smatch/smatch.c",
     "smatch/vector.c",
     "smatch/tree.c",
     "smatch/cat.c",
     "smatch/healpix.c"],
)
setup(
    name="smatch",
    packages=find_packages(),
    version="0.10.0",
    ext_modules=[ext],
    include_dirs=numpy.get_include(),
)
