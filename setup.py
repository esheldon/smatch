from setuptools import setup, Extension
import numpy


ext = Extension(
    "smatch._smatch",
    [
        "smatch/smatch.c",
        "smatch/vector.c",
        "smatch/tree.c",
        "smatch/cat.c",
        "smatch/healpix.c",
    ])

setup(
    name="smatch",
    packages=['smatch'],
    version="0.9.0.1",
    ext_modules=[ext],
    include_dirs=numpy.get_include(),
)
