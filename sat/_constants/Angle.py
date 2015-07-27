"""Necessary constants for tle.Angle class."""

from math import pi


# Angle modes.
DEG = "deg"
GRD = "grd"
RAD = "rad"
REV = "rev"

# Angle conversion constants.
DEG2RAD = pi/180
GRD2RAD = pi/200
RAD2RAD = 1.
REV2RAD = 2*pi

# Dictionary of constants and default values.
DEFAULT_ANGLE_VALUE = 0
INFINITY = float("inf")
PI = pi
