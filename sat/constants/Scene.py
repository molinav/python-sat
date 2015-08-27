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
DEFAULT_SCAN_STEP_LOOK_ANGLE =\
    0.053825

# Margin value when obtaining the border mask of a gridded indices array.
HPRT_BORDER_MASK_MARGIN =\
    200

# Threshold values for cloud masking of HRPT images. They are chosen to
# overestimate the presence of clouds, as cloud mask is used to remove non-valid
# key points from key-feature extraction.
HRPT_CLOUD_THRESHOLD_B5 =\
    500

# Default steps when georeferencing and HRPT image into geographic coordinates.
HRPT_LATITUDE_STEP =\
    -0.01
HRPT_LONGITUDE_STEP =\
    +0.01

# Constants value which determine the spacecraft's direction.
SPACECRAFT_ASCENDING =\
    +1
SPACECRAFT_DESCENDING =\
    -1
