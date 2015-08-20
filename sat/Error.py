"""Custom errors for classes from python-sat."""


class SatError(Exception):

    def __init__(self, message):
        """Constructor of a SatError instance."""
        self.message = message

    def __str__(self):
        """Generic text structure for printing a SatError."""

        title = self.__class__.__name__
        return "{}: {}".format(title, self.message)


class EphemerisError(SatError):
    pass


class FilepathError(SatError):
    pass


class OrbitError(SatError):
    pass


class SceneError(SatError):
    pass
