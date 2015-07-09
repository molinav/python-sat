#! /usr/bin/env/python

import re
from datetime import datetime
from datetime import timedelta
from math import radians
from six import text_type
from . TLESet import TLESet
from . TLEError import TLEError


def from_list(liststr):
    """Return a TLESet instance from a list of TLE string lines.

    Arguments:

        liststr : list

            list of three string lines containing the satellite title
            and the two-line element set

    Returns:

        obj : tle.TLESet

            TLESet instance containing all the parameters of an orbit

    """

    try:
        obj = TLESet()
        # Verify number of lines and characters per line.
        if not isinstance(liststr, list):
            raise TLEError(1001)
        if len(liststr) != obj.LINE_NUMBER:
            raise TLEError(1002)
        if max(not isinstance(item, str) for item in liststr):
            raise TLEError(1003)
        # Remove \n characters from string lines.
        liststr = [line.replace("\n", "") for line in liststr]
        if len(liststr[0]) > obj.TITLE_LENGTH:
            raise TLEError(1004)
        obj.satellite_name = text_type(liststr.pop(0).strip())
        if max(len(item) != obj.LINE_LENGTH for item in liststr):
            raise TLEError(1005)
        # Verify checksum values for both lines.
        try:
            checksum1 = [int(line[-1]) for line in liststr]
            checksum2 = [line.replace("-", "1") for line in liststr]
            checksum2 = [re.findall("[1-9]", line[:-1]) for line in checksum2]
            checksum2 = [sum(int(x) for x in line) % 10 for line in checksum2]
            if checksum1 != checksum2:
                raise ValueError
        except ValueError:
            err = TLEError(1006)
            err.__cause__ = None
            raise err
        # Verify start values for both lines.
        startnum = [line[0] for line in liststr]
        if startnum != ["1", "2"]:
            raise TLEError(1007)
        line1 = liststr[0]
        line2 = liststr[1]
        # Extract satellite number and classification.
        try:
            obj.satellite_number = int(line1[2:7])
            if obj.satellite_number != int(line2[2:7]):
                raise TLEError(1008)
        except ValueError:
            err = TLEError(1009)
            err.__cause__ = None
            raise err
        obj.satellite_classification = text_type(line1[7])
        # Extract satellite international designator (launch date and
        # launch item code).
        try:
            obj.international_designator = text_type(line1[9:17].strip())
            obj.launch_date = datetime.strptime(
                obj.international_designator[:-1], "%y%j").date()
            obj.launch_piece = text_type(obj.international_designator[-1])
        except ValueError:
            err = TLEError(1010)
            err.__cause__ = None
            raise err
        # Extract epoch year and fractional day.
        try:
            epoch_int = datetime.strptime(line1[18:23], "%y%j")
            epoch_dec = timedelta(days=float(line1[23:32]))
            obj.epoch = epoch_int + epoch_dec
        except ValueError:
            err = TLEError(1011)
            err.__cause__ = None
            raise err
        # Extract first time derivative of mean motion.
        try:
            obj.mean_motion_first_dif = 2 * float(line1[33:43])
        except ValueError:
            err = TLEError(1012)
            err.__cause__ = None
            raise err
        # Extract second time derivative of mean motion.
        try:
            sign = line1[44]
            mant = line1[45:50]
            expn = line1[50:52]
            num = float("{}.{}e{}".format(sign, mant, expn))
            obj.mean_motion_second_dif = 6 * num
        except ValueError:
            err = TLEError(1013)
            err.__cause__ = None
            raise err
        # Extract BSTAR aerodynamic drag term on the satellite.
        try:
            sign = line1[53]
            mant = line1[54:59]
            expn = line1[59:61]
            num = float("{}.{}e{}".format(sign, mant, expn))
            obj.bstar_drag = num
        except ValueError:
            err = TLEError(1014)
            err.__cause__ = None
            raise err
        # Extract ephemeris type (0 is always expected).
        try:
            obj.ephemeris_type = int(line1[62])
        except ValueError:
            err = TLEError(1015)
            err.__cause__ = None
            raise err
        # Extract element set number.
        try:
            obj.element_set_number = int(line1[64:68])
        except ValueError:
            err = TLEError(1016)
            err.__cause__ = None
            raise err
        # Extract inclination in degrees and convert to radians.
        try:
            obj.inclination = radians(float(line2[8:16]))
        except ValueError:
            err = TLEError(1017)
            err.__cause__ = None
            raise err
        # Extract right ascension of the ascending node in degrees and
        # convert to radians.
        try:
            obj.right_ascension_of_the_ascending_node = radians(
                float(line2[17:25]))
        except ValueError:
            err = TLEError(1018)
            err.__cause__ = None
            raise err
        # Extract orbital eccentricity.
        try:
            obj.eccentricity = float(".{}".format(line2[26:33]))
        except ValueError:
            err = TLEError(1019)
            err.__cause__ = None
            raise err
        # Extract argument of perigee in degrees and convert to radians.
        try:
            obj.argument_of_perigee = radians(float(line2[34:42]))
        except ValueError:
            err = TLEError(1020)
            err.__cause__ = None
            raise err
        # Extract mean anomaly in degrees and convert to radians.
        try:
            obj.mean_anomaly = radians(float(line2[43:51]))
        except ValueError:
            err = TLEError(1021)
            err.__cause__ = None
            raise err
        # Extract mean motion in revolutions per day.
        try:
            obj.mean_motion = float(line2[52:63])
        except ValueError:
            err = TLEError(1022)
            err.__cause__ = None
            raise err
        # Extract revolution number at epoch.
        try:
            obj.epoch_revolution_number = int(line2[63:68])
        except ValueError:
            err = TLEError(1023)
            err.__cause__ = None
            raise err
        # Return TLESet instance.
        return obj
    except TLEError as e:
        print(e)
        exit()
