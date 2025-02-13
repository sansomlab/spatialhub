import sys
import os
import re
import setuptools
from setuptools import setup, find_packages, Extension

from packaging.version import Version
if Version(setuptools.__version__) < Version('1.1'):
    print("Version detected:", Version(setuptools.__version__))
    raise ImportError(
        "spatialhub requires setuptools 1.1 higher")

########################################################################
########################################################################

IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################

# collect version
print(sys.path.insert(0, "spatialhub"))
import spatialhub.version as version
version = version.__version__

###############################################################
###############################################################

# Define dependencies
major, minor1, minor2, s, tmp = sys.version_info

if major < 3:
    raise SystemExit("""Requires Python 3 or later.""")

spatialhub_packages = find_packages()
spatialhub_package_dirs = {'spatialhub': 'spatialhub'}

##########################################################
##########################################################

# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='spatialhub',
    version=version,
    description='spatialhub: pipelines for spatial transcriptomics data pre-processing and analysis',
    author='Sansom lab',
    author_email='stephen.sansom@kennedy.ox.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='''spatialhub: pipelines for spatial transcriptomics data pre-processing and analysis''',
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="",
    # package contents
    packages=spatialhub_packages,
    package_dir=spatialhub_package_dirs,
    include_package_data=True,
    entry_points={
        "console_scripts": ["spatialhub = spatialhub.entry:main"]
    },
    # other options
    zip_safe=False,
    #test_suite="tests",
)