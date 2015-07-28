"""Custom errors for sat.Filepath class."""


class FilepathError(Exception):

    def __init__(self, message):
        """Constructor of an FilepathError instance."""
        self.message = message

    def __str__(self):
        """Generic text structure for printing a FilepathError."""

        title = self.__class__.__name__
        return "{}: {}".format(title, self.message)
