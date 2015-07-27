"""Private functions for nicely printing of TypeError messages."""


def str_format_types(types):
    return "({})".format(", ".join(t.__name__ for t in types))


def str_type_error_message(f, args, types, case):
    msg = "'{}' {} {}, but {} {}".format(
        f.__name__,
        ("accepts", "returns")[case],
        str_format_types(types),
        ("was given", "result is")[case],
        str_format_types(map(type, args)),
        )
    return msg
