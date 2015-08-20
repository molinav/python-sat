"""Example file where two Scene instances are created (day and night cases)
and geolocation coordinates in geographic WGS4 are computed."""

from __future__ import division
from __future__ import print_function
from six import text_type
from sat.Ephemeris import Ephemeris
from sat.Scene import Scene
import numpy as np
import os


# Define a conversion array from [rad, rad, m] into [deg, deg, km].
FACTOR = np.array([180/np.pi, 180/np.pi, 1/1000])
# Define output folder.
FOLDER = "./example3"

# Create output folder if it does not exist.
if not os.path.exists(FOLDER):
    os.mkdir(FOLDER)

# Open TLE file and extract its lines.
path_tle = "example3.tle"
with open(path_tle, "r") as file_tle:
    tle_list = [text_type(row.strip("\n")) for row in file_tle.readlines()]
ephem = Ephemeris.from_tle(*tle_list)

for case in ["day", "night"]:

    # Open MetOp-B Scene instance.
    path_img = text_type("example3_{}.l1b".format(case))
    scen = Scene(path_img, ephem)

    # Show time and TLE information.
    print("\n{}\n{}".format(path_img, 72*"-"))
    print("\nMetOp-B scene datetime information\n")
    print("  {:35} {}".format("Scene start datetime:    ", scen.start_datetime))
    print("  {:35} {}".format("Scene end datetime:      ", scen.end_datetime))
    print("  {:35} {}".format("Scene timedelta per scan:", scen.scan_timedelta))
    print(scen.orbit.ephemeris)

    # Extract original GCPs and save to GCP file.
    gcps = np.asarray([
        [item.GCPY, item.GCPX, item.GCPZ, item.GCPLine, item.GCPPixel]
        for item in scen.dataset.GetGCPs()])
    np.savetxt(os.path.join(FOLDER, "example3_{}_ref.gcp".format(case)),
               gcps, header="y,x,z,r,c", fmt="%.6f", delimiter=",")

    # Compute the appropriate coordinates for the Scene instance.
    # The spacing between GCP points is 40 pixels in x-direction (across-track)
    # and 23 pixels in y-direction (along-track).
    scen.compute(40, 23)
    trajectory = FACTOR * scen.orbit.position_geo
    coordtable = np.hstack([FACTOR * scen.position_geo, scen.position_pix])

    # Export trajectory and calculated GCPs grid into CSV files.
    # Header: (y, x, z)=(latitude, longitude, altitude), (r, c)=(row, column).
    np.savetxt(os.path.join(FOLDER, "example3_{}_trj.gcp".format(case)),
               trajectory, header="y,x,z", fmt="%.6f", delimiter=",")
    np.savetxt(os.path.join(FOLDER, "example3_{}_img.gcp".format(case)),
               coordtable, header="y,x,z,r,c", fmt="%.6f", delimiter=",")

    # Reformat table of coordinates to export into ENVI GCP file. Remember
    # that Python indexing starts in 0, while ENVI/IDL indexing starts in 1.
    coordtable_envi = coordtable[:, [1, 0, 4, 3]]
    coordtable_envi[:, 2:] += 1
    np.savetxt(os.path.join(FOLDER, "example3_{}_img.pts".format(case)),
               coordtable_envi,
               header="ENVI Rigorous Orthorectification GCP File",
               fmt="%.6f", delimiter="\t", comments=";")

    # Show end message.
    print("\nScene computation was run successfully\n")
