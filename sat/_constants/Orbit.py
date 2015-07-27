"""Necessary constants for tle.Orbit class."""

from six import text_type


# Limit values for Orbit properties.
LIMITS_ECCENTRICITY =\
    ["[", 0.0, 1.0, "]"]
LIMITS_ELEMENT_SET_NUMBER =\
    ["[", 0, 9999, "]"]
LIMITS_EPOCH_REVOLUTION_NUMBER =\
    ["[", 0, 99999, "]"]
LIMITS_MEAN_MOTION_FIRST_DIF =\
    ["[", -1., +1., "]"]
LIMITS_SATELLITE_NUMBER =\
    ["[", 0, 99999, "]"]

# Patterns used within the class Orbit.
PATTERN_LAUNCH_PIECE =\
    text_type("^([A-Z]{1,3})$")
PATTERN_SATELLITE_CLASSIFICATION =\
    text_type("^([A-Z])$")
PATTERN_SATELLITE_NAME =\
    text_type("^[^\t\n\r\f\v]{1,24}$")

# Patterns for a two-line element set.
PATTERN_TITLE = text_type(
    "^[^a-z\t\n\r\f\v]{24}$")
PATTERN_LINE1 = text_type(
    "^(1 {i5}{cap} {i5}{cod} {day} {dec} {ex} {ex} 0 {end}{chk})$".
    format(cap="([A-Z])",
           chk="(\d)",
           cod="([A-Z]  |[A-Z]{2} |[A-Z]{3})",
           day="(\d{2}([0-2]\d{2}|3[0-6]\d|36[0-5])\.\d{8})",
           dec="([ |-]\.\d{8})",
           end="(\d{4}| \d{3}|  \d{2}|   \d)",
           ex="([ |-]\d{5}[+|-]\d)",
           i5="(\d{5})",))
PATTERN_LINE2 = text_type(
    "^(2 {i5} {an1} {an2} {i7} {an2} {an2} {end}{chk})$".
    format(an1="((  \d|[ 1-2]\d{2}|3[0-5]\d).\d{4})",
           an2="((  \d| \d{2}|\d{3}).\d{4})",
           chk="(\d)",
           end="(( \d|\d{2})\.\d{8})(    \d|   \d{2}|  \d{3}| \d{4}|\d{5})",
           i5="(\d{5})",
           i7="(\d{7})"))
