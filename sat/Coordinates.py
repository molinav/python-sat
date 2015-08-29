from . constants.Coordinates import EARTH_SEMIMAJOR_AXIS
from . constants.Coordinates import EARTH_FLATTENING_FACTOR
from . constants.Coordinates import GEODETIC_COORDINATES_TOLERANCE
import numpy as np


class Coordinates(object):

    def __init__(self):
        pass

    @classmethod
    def from_ecf_to_geo(cls, obj):
        """Compute geodetic coordinates for a specific datetime.

        Parameters:

        obj
            (n, 3) array containing (x_ecef, y_ecef, z_ecef) in meters
        """

        # Call necessary constants.
        ae = EARTH_SEMIMAJOR_AXIS
        e2 = EARTH_FLATTENING_FACTOR
        # Call necessary properties.
        rx_s, ry_s, rz_s = [x[:, None] for x in obj.T]
        p_factor = np.sqrt(rx_s**2 + ry_s**2)

        def _estimate_n_factor(latitude):
            return ae / np.sqrt(1 - e2 * np.sin(latitude)**2)

        def _estimate_height(n_factor, latitude):
            return p_factor / np.cos(latitude) - n_factor

        def _estimate_latitude(n_factor, altitude):
            sin_lat = rz_s / (n_factor*(1-e2) + altitude)
            return np.arctan((rz_s + e2*n_factor*sin_lat)/p_factor)

        def _estimate_longitude():
            return np.arctan2(ry_s, rx_s)

        # Estimate recursively the values of altitude and latitude.
        n_factor_old = ae
        alt_old = 0
        lat_old = _estimate_latitude(n_factor_old, alt_old)
        dif = np.asarray([1])
        # Repeat recursive method until a tolerance is reached.
        while (dif > 0).all():
            n_factor_new = np.where(
                dif > 0, _estimate_n_factor(lat_old), n_factor_old)
            alt_new = np.where(
                dif > 0, _estimate_height(n_factor_new, lat_old), alt_old)
            lat_new = np.where(
                dif > 0, _estimate_latitude(n_factor_new, alt_new), lat_old)
            dif1 = np.abs(alt_new - alt_old)
            dif2 = np.abs(lat_new - lat_old)
            dif = dif1 + dif2 - GEODETIC_COORDINATES_TOLERANCE
            # Overwrite old temporary variables with new values.
            n_factor_old = n_factor_new
            alt_old = alt_new
            lat_old = lat_new
        # Compute longitude for a specific time.
        alt = alt_old
        lat = lat_old
        lon = _estimate_longitude()
        # Set geodetic satellite position property.
        return np.hstack([lat, lon, alt])

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
