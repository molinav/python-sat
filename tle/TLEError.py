#! /usr/bin/env/python


class TLEError(Exception):

    case = {
        1000: "undetermined error",
        1001: "input argument must be a list",
        1002: "input argument must contain three lines",
        1003: "input argument lines must be of type string",
        1004: "maximum title line length must be 24",
        1005: "maximum line length must be 69",
        1006: "at least one checksum is not valid",
        1007: "first value of a line must be the line number",
        1008: "satellite number do not match in both lines",
        1009: "invalid input for satellite number",
        1010: "invalid input for international designator",
        1011: "invalid input for epoch time",
        1012: "invalid input for first time derivative of mean motion",
        1013: "invalid input for second time derivative of mean motion",
        1014: "invalid input for BSTAR drag term",
        1015: "invalid input for ephemeris type",
        1016: "invalid input for element set number",
        1017: "invalid input for inclination angle",
        1018: "invalid input for right ascension of the ascending node",
        1019: "invalid input for eccentricity",
        1020: "invalid input for argument of perigee",
        1021: "invalid input for mean anomaly",
        1022: "invalid input for mean motion",
        1023: "invalid input for revolution number at epoch",
        1024: "input argument must be a file",
        }

    def __init__(self, code):
        """Constructor of a TLEError instance with a specific code."""

        try:
            self.code = code
            self.text = self.case[code]
        except KeyError:
            self.code = 1000
            self.text = self.case[code]

    def __str__(self):
        """Generic text structure for printing a TLEError."""

        title = self.__class__.__name__
        return "\n{}: {}".format(title, self.text)
