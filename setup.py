#! /usr/bin/env/python

from distutils.core import setup

__author__ = "Victor Molina"
__license__ = "GPLv3"
__version__ = "0.2dev"
__date__ = "2015/08/29"
__maintainer__ = "Victor Molina"
__email__ = "victor@latuv.uva.es"
__status__ = "Development"


setup(
    name="python-sat",
    version="0.1dev",
    packages=[
        "sat",
        "sat.constants",
        "sat.constans._matcher",
        "sat.decorators",
        "sat.filepath"],
    license="GNU General Public License Version 3",
    long_description=open("README.md").read(),
    install_requires=[
        "cv2",
        "gdal",
        "numpy",
        "scipy",
        "six",
        ],
    )
