import re
from six import text_type


def pattern(regex, getter=False, static=False):
    """Decorator that checks valid patterns of a string attribute.

    Parameters:

    regex
        regular expression as string
    getter
        flag that sets getter mode if True, otherwise setter mode is used
        (default False)
    static
        boolean flag that sets the decorator for static functions
        (default False)
    """

    def decorator(f):

        assert isinstance(regex, text_type)

        # Define inner function.
        def new_f(*args):
            # Assert characteristics of input arguments.
            flag = bool(not static)
            if getter:
                assert len(args) is 0 + flag
                val = f(*args)
            else:
                assert len(args) is 1 + flag
                val = args[flag]
            assert isinstance(val, text_type)
            # Check that the string variable has the expected length.
            if not re.match(regex, val):
                att = f.__name__.replace("_", " ")
                msg = "invalid {} '{}'".format(att, val)
                raise ValueError(msg)
            return val if getter else f(*args)

        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f

    # Return decorator.
    return decorator
