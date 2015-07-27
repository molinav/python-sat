from __future__ import print_function
import os.path
from six import text_type
from sat.Orbit import Orbit
from sat.Angle import Angle

expath = "example.tle"

with open(expath, "r") as exfile:

    tle_list1 = [text_type(exfile.readline().strip("\n")) for i in range(3)]
    
    print()
    for row in tle_list1:
        print(row)
    print()
    
    obj = Orbit.from_tle(*tle_list1)

    for attr in obj._properties:
        v = obj.__getattribute__(attr)
        if isinstance(v, Angle):
            v = v.deg
        print("{:40s}".format(attr), "{}".format(v))

    tle_list2 = obj.to_tle()
    
    flag = (tle_list1 == tle_list2)
    print("\nImported and exported TLE sets are equivalent: {}".format(flag))
