"""Necessary constants for sat.filepath.Filepath class and subclasses."""

from six import text_type


# Patterns for Filepath class.
PATTERN_FILEPATH_SATELLITE_NAME = text_type(
    "([-a-z]+)")
PATTERN_FILEPATH_PRODUCT = text_type(
    "([A-Z0-9\-]*)")

# Limits for FilepathNPP class.
LIMITS_FILEPATH_NPP_ORBIT_NUMBER =\
    ["[", 0, 99999, "]"]

# Patterns for FilepathNPP class.
PATTERN_FILEPATH_NPP_DATE = text_type(
    "((19|20)\d{2}(0\d|1[0-2])([0-2]\d|3[0-1]))")
PATTERN_FILEPATH_NPP_TIME = text_type(
    "(([0-1]\d|2[0-3])([0-5]\d)([0-5]\d)\d)")
PATTERN_FILEPATH_NPP_TIME_FULL = text_type(
    "(([0-1]\d|2[0-3])([0-5]\d)([0-5]\d)\d{6})")
PATTERN_FILEPATH_NPP_ORBIT_NUMBER = text_type(
    "(\d{5})")
PATTERN_FILEPATH_NPP_PROVIDER = text_type(
    "([a-z]{4})")
PATTERN_FILEPATH_NPP_STATUS = text_type(
    "([a-z]{3})")
PATTERN_FILEPATH_NPP = text_type(
    "^({pd}\_{sa}\_d{dt}\_t{t1}\_e{t1}\_b{ob}\_c{dt}{t2}\_{pr}\_{st}\.h5)$".
    format(dt=PATTERN_FILEPATH_NPP_DATE,
           ob=PATTERN_FILEPATH_NPP_ORBIT_NUMBER,
           pd=PATTERN_FILEPATH_PRODUCT,
           pr=PATTERN_FILEPATH_NPP_PROVIDER,
           sa=PATTERN_FILEPATH_SATELLITE_NAME,
           st=PATTERN_FILEPATH_NPP_STATUS,
           t1=PATTERN_FILEPATH_NPP_TIME,
           t2=PATTERN_FILEPATH_NPP_TIME_FULL,))
