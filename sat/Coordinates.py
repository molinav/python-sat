from . constants.Coordinates import EARTH_SEMIMAJOR_AXIS
from . constants.Coordinates import EARTH_FLATTENING_FACTOR
import numpy as np


class Coordinates(object):

    def __init__(self):
        pass

    @classmethod
    def from_geo_to_ecf(cls, obj):
        """Compute ECF coordinates using geodetic coordinates."""

        # Call necessary constants.
        ae = EARTH_SEMIMAJOR_AXIS
        e2 = EARTH_FLATTENING_FACTOR
        # Call necessary properties.
        lat, lon, alt = [x[:, None] for x in obj.T]

        # Compute ECF coordinates.
        n = ae / np.sqrt(1 - e2 * np.sin(lat)**2)
        rx_s = (n + alt) * np.cos(lat) * np.cos(lon)
        ry_s = (n + alt) * np.cos(lat) * np.sin(lon)
        rz_s = (n * (1 - e2) + alt) * np.sin(lat)
        # Set ECF satellite position property
        return np.hstack([rx_s, ry_s, rz_s])
