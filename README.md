# python-sat

The library `python-sat` is a simple Python interface to correct geometry of
AVHRR imagery using its satellite's orbital parameters from two-line element
sets (TLE) and estimating the clock drift and attitude position by feature
matching with reference datasets. It is released under the GNU
GPL v3.0.

Prerequisites
-------------

It is needed at least:

* Python 2.6+, 3.2+
* OpenCV 2.4.8+
* GDAL 1.11+
* NumPy 1.7+
* SciPy 0.15+
* Six for Python 2/3 compatibility

Installation
------------

Download the library as a compressed file and unzip it. Open a terminal in the
directory where the uncompressed folder is located and type:

```
pip install ./python-sat
```

To uninstall it, just type in a shell:

```
pip uninstall python-sat
```

Once the library is correctly installed, it is possible to check the example
file `example1.py`, which opens a TLE file and transcripts it into an Ephemeris
instance.

Reporting bugs
--------------

This library is currently in development and its state is still alpha. For
reporting bugs, feel free to open an issue at
http://github.com/molinav/python-sat/issues.
