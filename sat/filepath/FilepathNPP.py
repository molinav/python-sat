from .. constants.Filepath import LIMITS_FILEPATH_NPP_ORBIT_NUMBER
from .. constants.Filepath import PATTERN_FILEPATH_NPP_PROVIDER
from .. constants.Filepath import PATTERN_FILEPATH_NPP_STATUS
from .. constants.Filepath import PATTERN_FILEPATH_NPP
from .. decorators import accepts
from .. decorators import limits
from .. decorators import pattern
from .. decorators import returns
from . Filepath import Filepath
from datetime import datetime
from datetime import timedelta
from six import text_type
import os.path
import re


class FilepathNPP(Filepath):

    _properties = [
        # Parent properties.
        "satellite_name",
        "product",
        "start_datetime",
        "end_datetime",
        "processing_datetime",
        # New properties.
        "orbit_number",
        "provider",
        "status",
        ]

    def __init__(self, *args):
        """Constructor of a generic FilepathNPP instance.

        Parameters:

        args
            if args length is 0, an empty instance is created, if it
            is 1 an instance from string path is created, otherwise
            an error is raised
        """

        # Stop if no argument is provided.
        if len(args) is 0:
            super(FilepathNPP, self).__init__()
        # Parse a Suomi NPP path string if only one argument is provided.
        elif len(args) is 1:
            path = args[0]
            if isinstance(path, text_type):
                self.from_path(path)
            else:
                msg = "argument for FilepathNPP constructor must be a string"
                raise AttributeError(msg)
        # Raise an error otherwise.
        else:
            msg = "too many arguments for FilepathNPP constructor"
            raise AttributeError(msg)

    @property
    @returns(int)
    @limits(*LIMITS_FILEPATH_NPP_ORBIT_NUMBER, getter=True)
    def orbit_number(self):
        """orbit number"""
        return self._orbit_number

    @orbit_number.setter
    @accepts(int)
    @limits(*LIMITS_FILEPATH_NPP_ORBIT_NUMBER)
    def orbit_number(self, val):
        self._orbit_number = val

    @property
    @returns(text_type)
    @pattern(PATTERN_FILEPATH_NPP_PROVIDER, getter=True)
    def provider(self):
        """provider of the dataset"""
        return self._provider

    @provider.setter
    @accepts(text_type)
    @pattern(PATTERN_FILEPATH_NPP_PROVIDER)
    def provider(self, val):
        self._provider = val

    @property
    @returns(text_type)
    @pattern(PATTERN_FILEPATH_NPP_STATUS, getter=True)
    def status(self):
        """status of the satellite"""
        return self._status

    @status.setter
    @accepts(text_type)
    @pattern(PATTERN_FILEPATH_NPP_STATUS)
    def status(self, val):
        self._status = val

    @accepts(text_type)
    def from_path(self, path):
        """Fill FilepathNPP properties using an appropriate string path.


        Parameters:

        path
            filepath as a string variable.
        """

        # Create a temporary FilepathNPP instance.
        obj = FilepathNPP()
        path_name = os.path.basename(path)
        # Verify the path pattern.
        match0 = re.match(PATTERN_FILEPATH_NPP, path_name)
        if not match0:
            msg = "path string does not have a correct structure"
            raise AttributeError(msg)
        # Set dataset product name.
        obj.product = match0.group(1)
        # Set satellite name.
        obj.satellite_name = match0.group(2)
        # Set start and end datetimes.
        sdate = match0.group(3)
        stime = datetime.strptime(
            "".join([sdate, match0.group(4), "00000"]), "%Y%m%d%H%M%S%f")
        etime = datetime.strptime(
            "".join([sdate, match0.group(5), "00000"]), "%Y%m%d%H%M%S%f")
        if etime < stime:
            etime += timedelta(days=1)
        obj.start_datetime = stime
        obj.end_datetime = etime
        # Set orbit number.
        obj.orbit_number = int(match0.group(6))
        # Set processing datetime.
        pdate = match0.group(7)
        ptime = datetime.strptime(
            "".join([pdate, match0.group(8)]), "%Y%m%d%H%M%S%f")
        obj.processing_datetime = ptime
        # Set provider name.
        obj.provider = match0.group(9)
        # Set dataset status.
        obj.status = match0.group(10)
        # Copy attributes from temporary FilepathNPP instance into self.
        for item in obj._properties:
            item_name = "_{}".format(item)
            self.__setattr__(item_name, obj.__getattribute__(item_name))

    @returns(text_type)
    def to_path(self):
        """Return the appropriate string path for a FilepathNPP instance."""

        if not self:
            msg = "FilepathNPP instance is not complete"
            raise AttributeError(msg)

        out = text_type("{}_{}_d{}_t{}_e{}_b{:05d}_c{}_{}_{}.h5".
                        format(
                            self.product,
                            self.satellite_name,
                            self.start_datetime.strftime("%Y%m%d"),
                            self.start_datetime.strftime("%H%M%S%f")[:-5],
                            self.end_datetime.strftime("%H%M%S%f")[:-5],
                            self.orbit_number,
                            self.processing_datetime.strftime("%Y%m%d%H%M%S%f"),
                            self.provider,
                            self.status,))
        return out
