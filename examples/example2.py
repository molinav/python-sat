"""Example file where an Orbit instance is created and propagated."""

from __future__ import division
from __future__ import print_function
from datetime import timedelta
from six import text_type
from sat.Ephemeris import Ephemeris
from sat.Orbit import Orbit
import numpy as np
import os


# Define a conversion array from [rad, rad, m] into [deg, deg, km].
FACTOR = np.array([180/np.pi, 180/np.pi, 1/1000])
# Define output folder.
FOLDER = "./example2"

# Create output folder if it does not exist.
if not os.path.exists(FOLDER):
    os.mkdir(FOLDER)

# Open TLE file and create an Ephemeris instance.
path_tle = "example2.tle"
with open(path_tle, "r") as file_tle:
    tle_list = [text_type(row.strip("\n")) for row in file_tle.readlines()]
ephem = Ephemeris.from_tle(*tle_list)
print(ephem)

# Create Orbit instance with the previous Ephemeris object as input.
orb = Orbit(ephem)

# Extract epoch datetime and create a list of future datetimes.
sdate = orb.ephemeris.epoch_datetime
delta = timedelta(seconds=1)
dates = np.asarray([sdate + n*delta for n in range(50000)])[:, None]

# Propagate the orbit using the list of future datetimes.
orb.compute(dates)
coord = FACTOR * orb.position_geo

# Add an empty line everytime a flip in longitude takes place.
flips = np.hstack([[True], np.abs(coord[1:, 1]-coord[:-1, 1]) > 300])[:, None]
oarray = np.hstack([dates, coord, flips])
i = len(oarray) - 1
while True:
    if oarray[i, 4] is True:
        oarray = np.vstack([
            oarray[:i], np.array([5*[""]]), oarray[i:]])
    else:
        i -= 1
        if i is 0:
            break

# Export trajectory into CSV files.
# Header: (x, y, z)=(longitude, latitude, altitude).
np.savetxt(os.path.join(FOLDER, "example2_orbit.csv"), oarray[:, [2, 1, 3]],
           header="x,y,z", fmt="%s", delimiter=",", comments="")
