from . _str import str_type_error_message


def returns(types):
    """Decorator that checks the type of function's return value.

    Parameters:

    types
        the expected type of the decorated function's return value
    """

    def decorator(f):

        # Define inner function.
        def new_f(*args):
            # Assert characteristics of input arguments.
            assert isinstance(types, type)
            # Check that the output value has the expected type.
            val = f(*args)
            if not isinstance(val, types):
                msg = str_type_error_message(f, (val,), (types,), 1)
                raise TypeError(msg)
            return val

        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f

    # Return decorator.
    return decorator
