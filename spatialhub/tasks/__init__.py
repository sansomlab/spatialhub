'''
tasks.py
========

The :mod:`tasks` module contains helper functions for pipeline tasks.

Core components:

* `parameters`_
* `setup`_
* `zarr`_  [some SpatialData specific slot TBC]

Pipeline specific components:

* `cellxgene`_ [TBC]

'''


# import core submodules into top-level namespace

from spatialhub.tasks.setup import *
from spatialhub.tasks.parameters import *
#from spatialhub.tasks.zarr import *