"""Necessary constants for sat.Orbit class."""


# Conversion factor from days to seconds.
DAY2SEC =\
    86400

# Anomaly tolerance in rad.
ECCENTRIC_ANOMALY_TOLERANCE =\
    1e-9

# Coefficient for the second zonal term, related to the oblateness of the
# Earth, expressed in m^5 per s^2.
SECOND_ZONAL_TERM_J2 =\
    1.7555e25

# Conversion factor from mean solar time to mean sidereal time.
SOL2SID =\
    1.00273790935

# Standard gravitational parameter, that is, the product of the gravitational
# constant G and the Earth's mass M, expressed in m^3 per s^2.
STANDARD_GRAVITATIONAL_PARAMETER =\
    3.986004418e14
