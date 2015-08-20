"""Necessary constants for sat.Scene class."""

from six import text_type


# Set tolerance for recursive processes within sat.Scene class.
DATASET_COORDINATE_TOLERANCE =\
    1e-9

# GDAL dataset driver for AVHRR images.
DATASET_DRIVER = text_type(
    "NOAA Polar Orbiter Level 1b Data Set")

# Date pattern for AVHRR images opened with GDAL.
DATE_PATTERN = text_type(
    "^year: (\d{4}), day: (\d{1,3}), millisecond: (\d{1,8})$")

# Default scan step angle for AVHRR sensor in degrees.
DEFAULT_SCAN_STEP_ANGLE =\
    0.0539
