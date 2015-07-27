from . _str import str_type_error_message


def accepts(*types, **kwargs):
    """Decorator that checks the types of function's arguments.

    Parameters:

    types
        the expected types of the decorated function's input values
    kwargs
        keyword arguments only accepts boolean flag "static" that
        sets the decorator for static functions (default False)
    """

    # Verify the structure of keyword arguments.
    kwargs_flag = list(kwargs.keys())
    assert kwargs_flag == [] or kwargs_flag == ["static"]
    try:
        static = kwargs["static"]
    except KeyError:
        static = False

    def decorator(f):

        # Define inner function.
        def new_f(*args):
            # Assert characteristics of input arguments.
            flag = bool(not static)
            assert len(args) == len(types) + flag
            assert min(map(isinstance, types, (type for t in types)))
            # Check that the input values have the expected type.
            if not min(map(isinstance, args[flag:], types)):
                msg = str_type_error_message(f, args[flag:], types, 0)
                raise TypeError(msg)
            return f(*args)

        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f

    # Return decorator.
    return decorator
