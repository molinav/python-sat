# python-sat

The library `python-sat` is a simple Python interface to extract orbit
parameters from two-line element (TLE) sets. It is released under the GNU
GPL v3.0.

Prerequisites
-------------

It is needed at least:
* Python 2.6+, 3.2+
* The "six" package for Python 2/3 compatibility

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
file `example.py`, which opens a TLE file and transcripts it into an Orbit
instance.

Reporting bugs
--------------

This library is currently in development and its state is still alpha. For
reporting bugs, feel free to open an issue at
http://github.com/molinav/python-sat/issues.
