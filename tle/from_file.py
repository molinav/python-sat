#! /usr/bin/env/python

from io import TextIOWrapper
from . from_list import from_list
from . TLEError import TLEError


def from_file(filestr):
    """Return a TLESet instance from a text file containing TLE sets.

    Arguments:

        filestr : file

            text file containing two-line element sets

    Yields:

        obj : TLESet

            TLESet instance containing all the parameters of an orbit

    """

    # Verify that filestr is a text stream (string/unicode).
    try:
        if not isinstance(filestr, TextIOWrapper):
            try:
                if not isinstance(filestr, file):
                    raise TypeError
            except (NameError, TypeError):
                raise TLEError(1024)
        if filestr.mode != "r":
            message = "file cannot be read as text stream"
            raise OSError(message)
    except TLEError as e:
        print(e)
        exit()
    # Read all the lines within the file.
    liststr = filestr.readlines()
    if len(liststr) % 3 != 0:
        message = "number of file lines must be a multiple of 3"
        raise OSError(message)
    # Yield a TLESet instance for every group of three lines.
    num = len(liststr) // 3
    if num > 0:
        for i in range(0, num):
            obj = from_list(liststr[3*i:3*(i+1)])
            yield obj
