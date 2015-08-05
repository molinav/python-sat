"""Custom errors for sat.Ephemeris class."""


class EphemerisError(Exception):

    def __init__(self, message):
        """Constructor of an EphemerisError instance."""
        self.message = message

    def __str__(self):
        """Generic text structure for printing an EphemerisError."""

        title = self.__class__.__name__
        return "{}: {}".format(title, self.message)
