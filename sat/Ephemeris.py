from __future__ import division
from . _constants.Ephemeris import LIMITS_ECCENTRICITY
from . _constants.Ephemeris import LIMITS_ELEMENT_SET_NUMBER
from . _constants.Ephemeris import LIMITS_EPOCH_REVOLUTION_NUMBER
from . _constants.Ephemeris import LIMITS_MEAN_MOTION_FIRST_DIF
from . _constants.Ephemeris import LIMITS_SATELLITE_NUMBER
from . _constants.Ephemeris import PATTERN_LAUNCH_PIECE
from . _constants.Ephemeris import PATTERN_SATELLITE_CLASSIFICATION
from . _constants.Ephemeris import PATTERN_SATELLITE_NAME
from . _constants.Ephemeris import PATTERN_TITLE
from . _constants.Ephemeris import PATTERN_LINE1
from . _constants.Ephemeris import PATTERN_LINE2
from . _decorators import accepts
from . _decorators import limits
from . _decorators import pattern
from . _decorators import returns
from . _errors.Ephemeris import EphemerisError
from . Angle import AziAngle
from . Angle import ZenAngle
from datetime import datetime
from datetime import date
from datetime import timedelta
from six import text_type
import re


class Ephemeris(object):

    __slots__ = [
        "_argument_of_perigee",
        "_drag_term",
        "_eccentricity",
        "_element_set_number",
        "_epoch_datetime",
        "_epoch_revolution_number",
        "_inclination",
        "_launch_date",
        "_launch_piece",
        "_mean_anomaly",
        "_mean_motion",
        "_mean_motion_first_dif",
        "_mean_motion_second_dif",
        "_longitude_of_the_ascending_node",
        "_satellite_classification",
        "_satellite_name",
        "_satellite_number",
        ]

    def __init__(self):
        """Constructor of a generic Ephemeris instance."""

        for item in self.__slots__:
            self.__setattr__(item, None)

    def __bool__(self):
        """Return True if all the instance attributes are well-defined."""

        try:
            flag = min(self.__getattribute__(item) is not None
                       for item in self.__slots__)
        except (AssertionError, TypeError, ValueError):
            flag = False
        return flag

    def __repr__(self):
        """Fancy representation of an Ephemeris instance."""

        suplist = [
            [item[1:].replace("_", " ") + ":", self.__getattribute__(item)]
            for item in self.__slots__
            ]
        suplist = [
            [x, y if not isinstance(y, (AziAngle, ZenAngle)) else y.deg]
            for (x, y) in suplist
            ]

        text = (
            "\nEphemeris information\n"
            "\n  {:35s} {:.4f} deg"
            "\n  {:35s} {:.4e} per Earth radii"
            "\n  {:35s} {:.7f}"
            "\n  {:35s} {}"
            "\n  {:35s} {}"
            "\n  {:35s} {}"
            "\n  {:35s} {:.4f} deg"
            "\n  {:35s} {}"
            "\n  {:35s} {}"
            "\n  {:35s} {:.4f} deg"
            "\n  {:35s} {:.8f} rev per day"
            "\n  {:35s} {:.4e} rev per day^2"
            "\n  {:35s} {:.4e} rev per day^3"
            "\n  {:35s} {:.4f} deg"
            "\n  {:35s} {}"
            "\n  {:35s} {}"
            "\n  {:35s} {}").format(
                *[item for sublist in suplist for item in sublist])
        return text

    @staticmethod
    @accepts(text_type, text_type, text_type, static=True)
    def from_tle(title, line1, line2):
        """Return an Ephemeris instance from a two-line element set.

        Parameters:

        title
            two-line element set title
        line1
            first line of the two-line element set
        line2
            second line of the two-line element set
        """

        try:

            # Remove \n characters from string lines.
            obj = Ephemeris()
            lst = [x.replace("\n", "") for x in (title, line1, line2)]

            # Verify title, line1 and line2 with regular expressions.
            match0 = re.match(PATTERN_TITLE, lst[0])
            match1 = re.match(PATTERN_LINE1, lst[1])
            match2 = re.match(PATTERN_LINE2, lst[2])

            if not match0:
                msg = "invalid structure for TLE title"
                raise EphemerisError(msg)
            if not match1:
                msg = "invalid structure for TLE line 1"
                raise EphemerisError(msg)
            if not match2:
                msg = "invalid structure for TLE line 2"
                raise EphemerisError(msg)

            # Start transcription of title.
            tmp0 = match0.group(0).strip()
            obj.satellite_name = tmp0

            # Verify checksums.
            chk1 = [int(match1.group(11)), int(match2.group(9))]
            chk2 = [obj._calc_checksum(row[:-1]) for row in lst[1:]]
            if chk1 != chk2:
                msg = "checksum error, found {}, expected {}".format(chk1, chk2)
                raise EphemerisError(msg)

            # Verify and transcript satellite number.
            tmp1 = int(match1.group(1))
            tmp2 = int(match2.group(1))
            if tmp1 != tmp2:
                msg = "satellite number within the lines does not match"
                raise EphemerisError(msg)
            obj.satellite_number = tmp1

            # Start transcription of line1.
            # Transcript satellite classification.
            tmp1 = match1.group(2)
            obj.satellite_classification = tmp1

            # Transcript launch date.
            tmp1 = datetime.strptime(match1.group(3), "%y%j").date()
            obj.launch_date = tmp1

            # Transcript launch piece.
            tmp1 = match1.group(4).strip()
            obj.launch_piece = tmp1

            # Transcript epoch datetime.
            day1 = match1.group(5)
            day2 = match1.group(6)
            tmp1 = datetime.strptime(day1, "%y%j") + timedelta(days=float(day2))
            obj.epoch_datetime = tmp1

            # Transcript first derivative of mean motion.
            tmp1 = 2 * float(match1.group(7))
            obj.mean_motion_first_dif = tmp1

            # Transcript second derivative of mean motion.
            tmp1 = match1.group(8)
            tmp1 = 6 * float("{}.{}e{}".format(tmp1[0], tmp1[1:6], tmp1[6:]))
            obj.mean_motion_second_dif = tmp1

            # Transcript drag term.
            tmp1 = match1.group(9)
            tmp1 = float("{}.{}e{}".format(tmp1[0], tmp1[1:6], tmp1[6:]))
            obj.drag_term = tmp1

            # Transcript element set number.
            tmp1 = int(match1.group(10))
            obj.element_set_number = tmp1

            # Start transcription of line2.
            # Transcript inclination.
            tmp2 = ZenAngle(deg=float(match2.group(2)))
            obj.inclination = tmp2

            # Transcript longitude of the ascending node.
            tmp2 = AziAngle(deg=float(match2.group(3)))
            obj.longitude_of_the_ascending_node = tmp2

            # Transcript eccentricity.
            tmp2 = float(".{}".format(match2.group(4)))
            obj.eccentricity = tmp2

            # Transcript argument of perigee.
            tmp2 = AziAngle(deg=float(match2.group(5)))
            obj.argument_of_perigee = tmp2

            # Transcript mean anomaly.
            tmp2 = AziAngle(deg=float(match2.group(6)))
            obj.mean_anomaly = tmp2

            # Transcript mean motion.
            tmp2 = float(match2.group(7))
            obj.mean_motion = tmp2

            # Transcript revolution number at epoch.
            tmp2 = int(match2.group(8))
            obj._epoch_revolution_number = tmp2

            # Return the Ephemeris instance.
            return obj

        except EphemerisError as err:
            err.__cause__ = None
            print(err)

    def to_copy(self):
        """Return a deep copy of the Ephemeris instance."""

        obj = Ephemeris()
        for item in self.__slots__:
            obj.__setattr__(item, self.__getattribute__(item))
        return obj

    def to_tle(self):
        """Return a two-line element set as a list."""

        if not self:
            msg = "Ephemeris instance is not complete"
            raise AttributeError(msg)

        title = self.satellite_name.ljust(24)
        line1 = text_type(
            "1 {:5d}{:1s} {:5s}{:3s} {:5s}{:9s} {:1s}{:9s} {:8s} {:8s} 0 {:4d}".
            format(self.satellite_number,
                   self.satellite_classification,
                   self.launch_date.strftime("%y%j"),
                   self.launch_piece.ljust(3),
                   self.epoch_datetime.strftime("%y%j"),
                   "{:10.8f}".format(
                       timedelta(0,
                                 self.epoch_datetime.second,
                                 self.epoch_datetime.microsecond,
                                 0,
                                 self.epoch_datetime.minute,
                                 self.epoch_datetime.hour,
                                 0,
                                 ).total_seconds()/86400)[1:],
                   "-" if self.mean_motion_first_dif < 0 else " ",
                   "{:10.8f}".format(abs(self.mean_motion_first_dif)/2)[1:],
                   " 00000-0" if self.mean_motion_second_dif == 0 else
                   "{: =06d}{:+d}".
                   format(*[[int(round(float(x)*10000)), int(y)+1] for x, y in
                            ["{:.5e}".format(
                             self.mean_motion_second_dif/6).split("e")]][0]),
                   " 00000-0" if self.drag_term == 0 else
                   "{: =06d}{:+d}".
                   format(*[[int(round(float(x)*10000)), int(y)+1] for x, y in
                            ["{:.5e}".format(
                             self.drag_term).split("e")]][0]),
                   self.element_set_number))
        line2 = text_type(
            "2 {:5d} {:8.4f} {:8.4f} {:7s} {:8.4f} {:8.4f} {:11.8f}{:5d}".
            format(self.satellite_number,
                   self.inclination.deg,
                   self.longitude_of_the_ascending_node.deg,
                   "{:9.7f}".format(self.eccentricity)[2:],
                   self.argument_of_perigee.deg,
                   self.mean_anomaly.deg,
                   self.mean_motion,
                   self.epoch_revolution_number))
        line1 = "".join([line1, text_type(self._calc_checksum(line1))])
        line2 = "".join([line2, text_type(self._calc_checksum(line2))])
        return [title, line1, line2]

    @property
    @returns(AziAngle)
    def argument_of_perigee(self):
        """argument of perigee stored as an AziAngle instance"""
        return self._argument_of_perigee

    @argument_of_perigee.setter
    @accepts(AziAngle)
    def argument_of_perigee(self, val):
        self._argument_of_perigee = val

    @argument_of_perigee.deleter
    def argument_of_perigee(self):
        self._argument_of_perigee = None

    @property
    @returns(float)
    def drag_term(self):
        """BSTAR radiation pressure coefficient in per Earth radii"""
        return self._drag_term

    @drag_term.setter
    @accepts(float)
    def drag_term(self, val):
        self._drag_term = val

    @drag_term.deleter
    def drag_term(self):
        self._drag_term = None

    @property
    @returns(float)
    @limits(*LIMITS_ECCENTRICITY, getter=True)
    def eccentricity(self):
        """orbital eccentricity"""
        return self._eccentricity

    @eccentricity.setter
    @accepts(float)
    @limits(*LIMITS_ECCENTRICITY)
    def eccentricity(self, val):
        self._eccentricity = val

    @eccentricity.deleter
    def eccentricity(self):
        self._eccentricity = None

    @property
    @returns(int)
    @limits(*LIMITS_ELEMENT_SET_NUMBER, getter=True)
    def element_set_number(self):
        """element set number, incremented when a new two-line element
        set is generated for the satellite"""
        return self._element_set_number

    @element_set_number.setter
    @accepts(int)
    @limits(*LIMITS_ELEMENT_SET_NUMBER)
    def element_set_number(self, val):
        self._element_set_number = val

    @element_set_number.deleter
    def element_set_number(self):
        self._element_set_number = None

    @property
    @returns(datetime)
    def epoch_datetime(self):
        """datetime for the current epoch"""
        return self._epoch_datetime

    @epoch_datetime.setter
    @accepts(datetime)
    def epoch_datetime(self, val):
        self._epoch_datetime = val

    @epoch_datetime.deleter
    def epoch_datetime(self):
        self._epoch_datetime = None

    @property
    @returns(int)
    @limits(*LIMITS_EPOCH_REVOLUTION_NUMBER, getter=True)
    def epoch_revolution_number(self):
        """revolution number corresponding to the epoch"""
        return self._epoch_revolution_number

    @epoch_revolution_number.setter
    @accepts(int)
    @limits(*LIMITS_EPOCH_REVOLUTION_NUMBER)
    def epoch_revolution_number(self, val):
        self._epoch_revolution_number = val

    @epoch_revolution_number.deleter
    def epoch_revolution_number(self):
        self._epoch_revolution_number = None

    @property
    @returns(ZenAngle)
    def inclination(self):
        """orbital inclination as an ZenAngle instance"""
        return self._inclination

    @inclination.setter
    @accepts(ZenAngle)
    def inclination(self, val):
        self._inclination = val

    @inclination.deleter
    def inclination(self):
        self._inclination = None

    @property
    @returns(date)
    def launch_date(self):
        """date of launch of the satellite"""
        return self._launch_date

    @launch_date.setter
    @accepts(date)
    def launch_date(self, val):
        self._launch_date = val

    @launch_date.deleter
    def launch_date(self):
        self._launch_date = None

    @property
    @returns(text_type)
    @pattern(PATTERN_LAUNCH_PIECE, getter=True)
    def launch_piece(self):
        """piece of launch which corresponds to the satellite"""
        return self._launch_piece

    @launch_piece.setter
    @accepts(text_type)
    @pattern(PATTERN_LAUNCH_PIECE)
    def launch_piece(self, val):
        self._launch_piece = val

    @launch_piece.deleter
    def launch_piece(self):
        self._launch_piece = None

    @property
    @returns(AziAngle)
    def longitude_of_the_ascending_node(self):
        """longitude of the ascending node stored as an AziAngle instance"""
        return self._longitude_of_the_ascending_node

    @longitude_of_the_ascending_node.setter
    @accepts(AziAngle)
    def longitude_of_the_ascending_node(self, val):
        self._longitude_of_the_ascending_node = val

    @longitude_of_the_ascending_node.deleter
    def longitude_of_the_ascending_node(self):
        self._longitude_of_the_ascending_node = None

    @property
    @returns(AziAngle)
    def mean_anomaly(self):
        """mean anomaly of the orbit stored as an AziAngle instance"""
        return self._mean_anomaly

    @mean_anomaly.setter
    @accepts(AziAngle)
    def mean_anomaly(self, val):
        self._mean_anomaly = val

    @mean_anomaly.deleter
    def mean_anomaly(self):
        self._mean_anomaly = None

    @property
    @returns(float)
    def mean_motion(self):
        """mean motion (time-average angular velocity) of the satellite
        over an orbit in revs per day"""
        return self._mean_motion

    @mean_motion.setter
    @accepts(float)
    def mean_motion(self, val):
        self._mean_motion = val

    @mean_motion.deleter
    def mean_motion(self):
        self._mean_motion = None

    @property
    @returns(float)
    @limits(*LIMITS_MEAN_MOTION_FIRST_DIF, getter=True)
    def mean_motion_first_dif(self):
        """first time derivative of the mean motion in revs per day^2"""
        return self._mean_motion_first_dif

    @mean_motion_first_dif.setter
    @accepts(float)
    @limits(*LIMITS_MEAN_MOTION_FIRST_DIF)
    def mean_motion_first_dif(self, val):
        self._mean_motion_first_dif = val

    @mean_motion_first_dif.deleter
    def mean_motion_first_dif(self):
        self._mean_motion_first_dif = None

    @property
    @returns(float)
    def mean_motion_second_dif(self):
        """second time derivative of the mean motion in revs per day^3"""
        return self._mean_motion_second_dif

    @mean_motion_second_dif.setter
    @accepts(float)
    def mean_motion_second_dif(self, val):
        self._mean_motion_second_dif = val

    @mean_motion_second_dif.deleter
    def mean_motion_second_dif(self):
        self._mean_motion_second_dif = None

    @property
    @returns(text_type)
    @pattern(PATTERN_SATELLITE_CLASSIFICATION, getter=True)
    def satellite_classification(self):
        """classification of the satellite (U = unclassified)"""
        return self._satellite_classification

    @satellite_classification.setter
    @accepts(text_type)
    @pattern(PATTERN_SATELLITE_CLASSIFICATION)
    def satellite_classification(self, val):
        self._satellite_classification = val

    @satellite_classification.deleter
    def satellite_classification(self):
        self._satellite_classification = None

    @property
    @returns(text_type)
    @pattern(PATTERN_SATELLITE_NAME, getter=True)
    def satellite_name(self):
        """name of the satellite"""
        return self._satellite_name

    @satellite_name.setter
    @accepts(text_type)
    @pattern(PATTERN_SATELLITE_NAME)
    def satellite_name(self, val):
        self._satellite_name = val

    @satellite_name.deleter
    def satellite_name(self):
        self._satellite_name = None

    @property
    @returns(int)
    @limits(*LIMITS_SATELLITE_NUMBER, getter=True)
    def satellite_number(self):
        """number code which identifies the satellite"""
        return self._satellite_number

    @satellite_number.setter
    @accepts(int)
    @limits(*LIMITS_SATELLITE_NUMBER)
    def satellite_number(self, val):
        self._satellite_number = val

    @satellite_number.deleter
    def satellite_number(self):
        self._satellite_number = None

    @accepts(text_type)
    def _calc_checksum(self, line):
        """Calculate the checksum mod 10 of a TLE line."""

        nums = re.findall("\d", line.replace("-", "1"))
        return sum(int(x) for x in nums) % 10 if nums else 0
