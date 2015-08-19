"""Custom errors for sat.Scene class."""


class SceneError(Exception):

    def __init__(self, message):
        """Constructor of a SceneError instance."""
        self.message = message

    def __str__(self):
        """Generic text structure for printing a SceneError."""

        title = self.__class__.__name__
        return "{}: {}".format(title, self.message)
