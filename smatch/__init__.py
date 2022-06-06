# flake8: noqa
from . import smatch
from .smatch import (
    match,
    match_self,
    Catalog,
    read_matches,
    match_dtype,
)
from .matcher import (
    Matcher,
    sphdist,
)

__version__ = '0.10.0'
