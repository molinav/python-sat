from six import text_type


def maximum_text_length(max_length, getter=False, static=False):
    """Decorator that checks the maximum length of a string value.

    Parameters:

    max_length
        maximum value for string length
    getter
        flag that sets getter mode if True, otherwise setter mode is used
        (default False)
    static
        boolean flag that sets the decorator for static functions
        (default False)
    """

    def decorator(f):

        assert isinstance(max_length, int)
        assert max_length > 0

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
            if len(val) > max_length:
                msg = "text string exceeds the maximum length ({})"\
                    .format(max_length)
                raise ValueError(msg)
            return val if getter else f(*args)

        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f

    # Return decorator.
    return decorator
