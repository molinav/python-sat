"""Necessary constants for sat.Coordinates class."""


# Earth's constants related to WGS84 ellipsoid (SI units).
EARTH_FLATTENING_VALUE =\
    1/298.257223563
EARTH_FLATTENING_FACTOR =\
    EARTH_FLATTENING_VALUE * (2 - EARTH_FLATTENING_VALUE)
EARTH_SEMIMAJOR_AXIS =\
    6378137
