#! /usr/bin/env/python

from datetime import date
from datetime import datetime
from math import pi
from six import text_type
from . decorators import accepts
from . decorators import maximum_text_length
from . decorators import number_limits


class TLESet(object):

    _properties = [
        "argument_of_perigee",
        "bstar_drag",
        "eccentricity",
        "element_set_number",
        "ephemeris_type",
        "epoch",
        "epoch_revolution_number",
        "inclination",
        "international_designator",
        "launch_date",
        "launch_piece",
        "mean_anomaly",
        "mean_motion",
        "mean_motion_first_dif",
        "mean_motion_second_dif",
        "right_ascension_of_the_ascending_node",
        "satellite_classification",
        "satellite_name",
        "satellite_number",
        ]

    LINE_NUMBER = 3
    LINE_LENGTH = 69
    TITLE_LENGTH = 24

    def __init__(self):
        """Constructor of a generic TLESet instance."""

        for item in self._properties:
            self.__setattr__("_{}".format(item), None)

    @property
    def argument_of_perigee(self):
        """argument of perigee in radians"""
        return self._argument_of_perigee

    @argument_of_perigee.setter
    @accepts(float)
    @number_limits(0, 2*pi)
    def argument_of_perigee(self, value):
        self._argument_of_perigee = value

    @argument_of_perigee.deleter
    def argument_of_perigee(self):
        self._argument_of_perigee = None

    @property
    def bstar_drag(self):
        """BSTAR radiation pressure coefficient in (Earth radii)^(-1)"""
        return self._bstar_drag

    @bstar_drag.setter
    @accepts(float)
    def bstar_drag(self, value):
        self._bstar_drag = value

    @bstar_drag.deleter
    def bstar_drag(self):
        self._bstar_drag = None

    @property
    def eccentricity(self):
        """orbital eccentricity"""
        return self._eccentricity

    @eccentricity.setter
    @accepts(float)
    @number_limits(0, 1)
    def eccentricity(self, value):
        self._eccentricity = value

    @eccentricity.deleter
    def eccentricity(self):
        self._eccentricity = None

    @property
    def element_set_number(self):
        """element set number, incremented when a new two-line element set
        is generated for the satellite"""
        return self._element_set_number

    @element_set_number.setter
    @accepts(int)
    @number_limits(0, 999)
    def element_set_number(self, value):
        self._element_set_number = value

    @element_set_number.deleter
    def element_set_number(self):
        self._element_set_number = None

    @property
    def ephemeris_type(self):
        """ephemeris type (actually 0 is always stored)"""
        return self._ephemeris_type

    @ephemeris_type.setter
    @accepts(int)
    @number_limits(0, 9)
    def ephemeris_type(self, value):
        self._ephemeris_type = value

    @ephemeris_type.deleter
    def ephemeris_type(self):
        self._ephemeris_type = None

    @property
    def epoch(self):
        """datetime for the current epoch time"""
        return self._epoch

    @epoch.setter
    @accepts(datetime)
    def epoch(self, value):
        self._epoch = value

    @epoch.deleter
    def epoch(self):
        self._epoch = None

    @property
    def epoch_revolution_number(self):
        """revolution number corresponding to the epoch time"""
        return self._epoch_revolution_number

    @epoch_revolution_number.setter
    @accepts(int)
    @number_limits(0, 99999)
    def epoch_revolution_number(self, value):
        self._epoch_revolution_number = value

    @epoch_revolution_number.deleter
    def epoch_revolution_number(self):
        self._epoch_revolution_number = None

    @property
    def inclination(self):
        """orbital inclination in radians"""
        return self._inclination

    @inclination.setter
    @accepts(float)
    @number_limits(0, pi)
    def inclination(self, value):
        self._inclination = value

    @inclination.deleter
    def inclination(self):
        self._inclination = None

    @property
    def international_designator(self):
        """international code which identifies the satellite"""
        return self._international_designator

    @international_designator.setter
    @accepts(text_type)
    @maximum_text_length(8)
    def international_designator(self, value):
        self._international_designator = value

    @international_designator.deleter
    def international_designator(self):
        self._international_designator = None

    @property
    def launch_date(self):
        """date of launch of the satellite"""
        return self._launch_date

    @launch_date.setter
    @accepts(date)
    def launch_date(self, value):
        self._launch_date = value

    @launch_date.deleter
    def launch_date(self):
        self._launch_date = None

    @property
    def launch_piece(self):
        """piece of launch which corresponds to the satellite"""
        return self._launch_piece

    @launch_piece.setter
    @accepts(text_type)
    @maximum_text_length(3)
    def launch_piece(self, value):
        self._launch_piece = value

    @launch_piece.deleter
    def launch_piece(self):
        self._launch_piece = None

    @property
    def mean_anomaly(self):
        """mean anomaly of the Kepler orbit in radians"""
        return self._mean_anomaly

    @mean_anomaly.setter
    @accepts(float)
    @number_limits(0, 2*pi)
    def mean_anomaly(self, value):
        self._mean_anomaly = value

    @mean_anomaly.deleter
    def mean_anomaly(self):
        self._mean_anomaly = None

    @property
    def mean_motion(self):
        """mean motion (time-average angular velocity) of the satellite
        over an orbit in revs per day"""
        return self._mean_motion

    @mean_motion.setter
    @accepts(float)
    def mean_motion(self, value):
        self._mean_motion = value

    @mean_motion.deleter
    def mean_motion(self):
        self._mean_motion = None

    @property
    def mean_motion_first_dif(self):
        """first time derivative of the mean motion in revs per day^2"""
        return self._mean_motion_first_dif

    @mean_motion_first_dif.setter
    @accepts(float)
    def mean_motion_first_dif(self, value):
        self._mean_motion_first_dif = value

    @mean_motion_first_dif.deleter
    def mean_motion_first_dif(self):
        self._mean_motion_first_dif = None

    @property
    def mean_motion_second_dif(self):
        """second time derivative of the mean motion in revs per day^3"""
        return self._mean_motion_second_dif

    @mean_motion_second_dif.setter
    @accepts(float)
    def mean_motion_second_dif(self, value):
        self._mean_motion_second_dif = value

    @mean_motion_second_dif.deleter
    def mean_motion_second_dif(self):
        self._mean_motion_second_dif = None

    @property
    def right_ascension_of_the_ascending_node(self):
        """right ascension of the ascending node in radians"""
        return self._right_ascension_of_the_ascending_node

    @right_ascension_of_the_ascending_node.setter
    @accepts(float)
    @number_limits(0, 2*pi)
    def right_ascension_of_the_ascending_node(self, value):
        self._right_ascension_of_the_ascending_node = value

    @right_ascension_of_the_ascending_node.deleter
    def right_ascension_of_the_ascending_node(self):
        self._right_ascension_of_the_ascending_node = None

    @property
    def satellite_classification(self):
        """classification of the satellite (U = unclassified)"""
        return self._satellite_classification

    @satellite_classification.setter
    @accepts(text_type)
    @maximum_text_length(1)
    def satellite_classification(self, value):
        self._satellite_classification = value

    @satellite_classification.deleter
    def satellite_classification(self):
        self._satellite_classification = None

    @property
    def satellite_name(self):
        """name of the satellite"""
        return self._satellite_name

    @satellite_name.setter
    @accepts(text_type)
    @maximum_text_length(24)
    def satellite_name(self, value):
        self._satellite_name = value

    @satellite_name.deleter
    def satellite_name(self):
        self._satellite_name = None

    @property
    def satellite_number(self):
        """number code which identifies the satellite"""
        return self._satellite_number

    @satellite_number.setter
    @accepts(int)
    @number_limits(0, 99999)
    def satellite_number(self, value):
        self._satellite_number = value

    @satellite_number.deleter
    def satellite_number(self):
        self._satellite_number = None
