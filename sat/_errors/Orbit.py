"""Custom errors for sat.Orbit class."""


class OrbitError(Exception):

    def __init__(self, message):
        """Constructor of an OrbitError instance."""
        self.message = message

    def __str__(self):
        """Generic text structure for printing an OrbitError."""

        title = self.__class__.__name__
        return "{}: {}".format(title, self.message)
