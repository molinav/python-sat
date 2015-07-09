#! /usr/bin/env/python

import tle


# Example 1.
input_list = [
    "METOP-B      \n",
    "1 38771U 12049A   15060.09038281  .00000104  00000-0  67552-4 0  9999\n",
    "2 38771  98.7116 121.5278 0001914 123.9555 292.2003 14.21473887127060\n",
]

sat1 = tle.from_list(input_list)

print("\nTest with tle.from_list function")
print("{}\n".format(sat1))

for k, v in sorted(sat1.__dict__.items()):
    print("{:40} {}".format(k[1:], v))

# Example 2.
path = "example.tle"
with open(path, "r") as input_file:
    sat2 = tuple(tle.from_file(input_file))[0]

print("\nTest with tle.from_file function")
print("{}\n".format(sat2))

for k, v in sorted(sat2.__dict__.items()):
    print("{:40} {}".format(k[1:], v))
