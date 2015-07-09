#! /usr/bin/env/python


def _info(f, expected, actual, flag):
    """Return nicely formatted error/warning message."""

    def format_types(types):
        return ", ".join([str(t).split("'")[1] for t in types])

    message = "'{}' {} ({}), but {} ({})".format(
        f.__name__,
        ("accepts", "returns")[flag],
        format_types(expected),
        ("was given", "result is")[flag],
        format_types(actual),
        )
    return message

def accepts(*types):
    """Decorator that checks types of function's arguments.

    Parameters:

    types
        the expected types of the inputs of the decorated function
    """

    def decorator(f):
        # Define inner function.
        def new_f(self, *args):
            assert len(args) == len(types)
            arg_types = tuple(map(type, args))
            if arg_types != types:
                message = _info(f, types, arg_types, 0)
                raise TypeError(message)
            return f(self, *args)
        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f
    # Return decorator.
    return decorator

def returns(return_type):
    """Decorator that checks types of function's return value.

    Parameters:

    return_type
        the expected type of the decorated function's return value
    """

    def decorator(f):
        # Define inner function.
        def new_f(self, *args):
            result = f(self, *args)
            result_type = type(result)
            if result_type != return_type:
                message = _info(f, (return_type,), (result_type,), 1)
                raise TypeError(message)
            return result
        new_f.__name__ = f.__name__
        # Return inner function.
        return new_f
    # Return decorator.
    return decorator

def maximum_text_length(max_length):
    """Decorator that checks the length of a string property.

    Parameters:

    max_length
        the maximum value for text length
    """

    def decorator(f):
        # Define inner function.
        def new_f(self, text):
            if len(text) > max_length:
                message = "maximum text length for '{}' is {}".format(
                    f.__name__,
                    max_length,
                    )
                raise ValueError(message)
            # Return inner function.
            return f(self, text)
        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f
    # Return decorator.
    return decorator

def number_limits(min_value, max_value):
    """Decorator that checks the belonging of a property to the interval
    defined by [min_value, max_value].

    Parameters:

    min_value
        the minimum accepted value (None means there is no minimum limit)
    max_value
        the maximum accepted value (None means there is no maximum limit)
    """

    def decorator(f):
        # Define inner function.
        def new_f(self, value):
            if min_value is not None:
                if value < min_value:
                    message = "minimum value for '{}' is {}".format(
                        f.__name__,
                        min_value,
                        )
                    raise ValueError(message)
            if max_value is not None:
                if value > max_value:
                    message = "maximum value for '{}' is {}".format(
                        f.__name__,
                        max_value,
                        )
                    raise ValueError(message)
            # Return inner function.
            return f(self, value)
        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f
    # Return decorator.
    return decorator
