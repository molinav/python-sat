from numbers import Real


def limits(left, minv, maxv, right, getter=False, static=False):
    """Getter decorator that checks the value is within the valid range.

    Parameters:

    left
        type of lower limit ('(' for open, '[' for closed)
    minv
        the lower limit value
    maxv
        the upper limit value
    right
        type of upper limit (')' for open, ']' for closed)
    getter
        flag that sets getter mode if True, otherwise setter mode is used
        (default False)
    static
        boolean flag that sets the decorator for static functions
        (default False)
    """

    def decorator(f):

        assert left in "([" and right in ")]"
        assert min(map(isinstance, (minv, maxv), (Real, Real)))
        assert minv < maxv

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
            assert isinstance(val, Real)
            # Check that input value is within the expected range.
            min_flag = val <= minv if left is "(" else val < minv
            max_flag = val >= maxv if left is ")" else val > maxv
            if min_flag or max_flag:
                msg = "value '{}' is not within range [{:.4g}, {:.4g}]"\
                    .format(val, minv, maxv)
                raise TypeError(msg)
            return val if getter else f(*args)

        # Return inner function.
        new_f.__name__ = f.__name__
        return new_f

    # Return decorator.
    return decorator
