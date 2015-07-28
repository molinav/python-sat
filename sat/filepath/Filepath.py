from .. _constants.Filepath import PATTERN_FILEPATH_SATELLITE_NAME
from .. _constants.Filepath import PATTERN_FILEPATH_PRODUCT
from .. _decorators import accepts
from .. _decorators import pattern
from .. _decorators import returns
from .. _errors.Filepath import FilepathError
from datetime import datetime
from six import text_type


class Filepath(object):

    _properties = [
        "satellite_name",
        "product",
        "start_datetime",
        "end_datetime",
        "processing_datetime",
        ]

    def __init__(self):
        """Constructor of a generic Filepath instance."""

        for item in self._properties:
            self.__setattr__("_{}".format(item), None)

    def __bool__(self):
        """Return True if all the instance attributes are well-defined."""

        try:
            flag = min(self.__getattribute__("_{}".format(item)) is not None
                       for item in self._properties)
        except (AssertionError, TypeError, ValueError):
            flag = False
        return flag

    @property
    @returns(text_type)
    @pattern(PATTERN_FILEPATH_PRODUCT, getter=True)
    def product(self):
        """name of the satellite product"""
        return self._product

    @product.setter
    @accepts(text_type)
    @pattern(PATTERN_FILEPATH_PRODUCT)
    def product(self, val):
        self._product = val

    @property
    @returns(text_type)
    @pattern(PATTERN_FILEPATH_SATELLITE_NAME, getter=True)
    def satellite_name(self):
        """name of the satellite"""
        return self._satellite_name

    @satellite_name.setter
    @accepts(text_type)
    @pattern(PATTERN_FILEPATH_SATELLITE_NAME)
    def satellite_name(self, val):
        self._satellite_name = val

    @property
    @returns(datetime)
    def start_datetime(self):
        """datetime for the first dataset scan line received"""

        try:
            if self._end_datetime < self._start_datetime:
                msg = "start datetime must be before end datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        try:
            if self._processing_datetime < self._start_datetime:
                msg = "start datetime must be before processing datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        return self._start_datetime

    @start_datetime.setter
    @accepts(datetime)
    def start_datetime(self, val):

        try:
            if self._end_datetime < val:
                msg = "start datetime must be before end datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        try:
            if self._processing_datetime < val:
                msg = "start datetime must be before processing datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        self._start_datetime = val

    @property
    @returns(datetime)
    def end_datetime(self):
        """datetime for the last dataset scan line received"""

        try:
            if self._start_datetime > self._end_datetime:
                msg = "end datetime must be after start datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        try:
            if self._processing_datetime < self._end_datetime:
                msg = "end datetime must be before processing datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        return self._end_datetime

    @end_datetime.setter
    @accepts(datetime)
    def end_datetime(self, val):

        try:
            if self._start_datetime > val:
                msg = "end datetime must be after start datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        try:
            if self._processing_datetime < val:
                msg = "end datetime must be before processing datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        self._end_datetime = val

    @property
    @returns(datetime)
    def processing_datetime(self):
        """datetime in which the dataset was processed"""

        try:
            if self._start_datetime > self._processing_datetime:
                msg = "processing datetime must be after start datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        try:
            if self._end_datetime > self._processing_datetime:
                msg = "processing datetime must be after end datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        return self._processing_datetime

    @processing_datetime.setter
    @accepts(datetime)
    def processing_datetime(self, val):

        try:
            if self._start_datetime > val:
                msg = "processing datetime must be after start datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        try:
            if self._end_datetime > val:
                msg = "processing datetime must be after end datetime"
                raise ValueError(msg)
        except TypeError:
            pass

        self._processing_datetime = val

    @accepts(text_type)
    def from_path(self, path):
        """Fill Filepath properties using an appropriate string path."""
    pass

    @returns(text_type)
    def to_path(self):
        """Return the appropriate string path for a Filepath instance."""
        pass

    @returns(bool)
    def is_within(self, other):
        """Return True if a Filepath instance is temporally within
        another Filepath instance."""

        if not isinstance(other, type(self)):
            msg = "argument must be of type {}".format(self.__class__.__name__)
            raise FilepathError(msg)
        if not self or not other:
            msg = "Path instances are not complete"
            raise FilepathError(msg)

        if self.start_datetime < other.start_datetime:
            return False
        elif self.end_datetime > other.end_datetime:
            return False
        else:
            return True
