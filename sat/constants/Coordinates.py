"""Necessary constants for sat.Coordinates class."""


# Earth's angular speed in radians per second.
EARTH_ANGULAR_SPEED =\
    7.2921150e-5

# Earth's constants related to WGS84 ellipsoid (SI units).
EARTH_FLATTENING_VALUE =\
    1/298.257223563
EARTH_FLATTENING_FACTOR =\
    EARTH_FLATTENING_VALUE * (2 - EARTH_FLATTENING_VALUE)
EARTH_SEMIMAJOR_AXIS =\
    6378137

# Coordinate tolerance for geodetic calculations.
GEODETIC_COORDINATES_TOLERANCE =\
    1e-9
