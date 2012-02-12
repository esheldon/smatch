import distutils
from distutils.core import setup, Extension, Command
import os
import numpy


ext=Extension("smatch._smatch", 
              ["smatch/smatch.c",
               "smatch/stack.c",
               "smatch/tree.c",
               "smatch/healpix.c"])
setup(name="smatch", 
      packages=['smatch'],
      version="0.1",
      ext_modules=[ext],
      include_dirs=numpy.get_include())




