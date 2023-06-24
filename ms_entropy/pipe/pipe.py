import pandas as pd
import functools

__all__ = ["Pipe", "pipe_function", "filter", "apply", "assign_and_apply"]


class Pipe(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __rrshift__(self, other):
        return Pipe(other)

    def __rshift__(self, other):
        func, args, kwargs = other
        return Pipe(func(self, *args, **kwargs))


def pipe_function(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func, args, kwargs

    return wrapper


@pipe_function
def filter(df, func, *args, **kwargs):
    f = lambda x: func(x, *args, **kwargs)
    return df[df.apply(f, axis=1)]


@pipe_function
def apply(df, func, *args, **kwargs):
    f = lambda x: func(x, *args, **kwargs)
    return df.apply(f, axis=1)


@pipe_function
def assign_and_apply(df, col, func, *args, **kwargs):
    f = lambda x: func(x, *args, **kwargs)
    return df.assign(**{col: df.apply(f, axis=1)})
