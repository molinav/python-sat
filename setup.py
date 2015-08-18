#! /usr/bin/env/python

from distutils.core import setup

__author__ = "Victor Molina"
__license__ = "GPLv3"
__version__ = "0.1dev"
__date__ = "2015/07/09"
__maintainer__ = "Victor Molina"
__email__ = "victor@latuv.uva.es"
__status__ = "Development"


setup(
    name="python-sat",
    version="0.1dev",
    packages=[
        "sat",
        "sat._constants",
        "sat._decorators",
        "sat._errors",
        "sat.filepath"],
    license="GNU General Public License Version 3",
    long_description=open("README.md").read(),
    install_requires=[
        "numpy",
        "six",
        ],
    )
