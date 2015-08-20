"""Example file where an Ephemeris instance is created using a TLE file."""

from __future__ import print_function
from six import text_type
from sat.Ephemeris import Ephemeris


path_tle = "example1.tle"

with open(path_tle, "r") as file_tle:

    # Read TLE text, create Ephemeris instance and export again to TLE text.
    tle_list1 = [text_type(row.strip("\n")) for row in file_tle.readlines()]
    obj = Ephemeris.from_tle(*tle_list1)
    tle_list2 = obj.to_tle()

    # Check that the original TLE and the new one are equivalent.
    flag = (tle_list1 == tle_list2)

    # Print TLE information.
    print("\n{}\n{}\n{}".format(*tle_list1))
    print(obj)
    print("\nImported and exported TLE sets are equivalent: {}".format(flag))
