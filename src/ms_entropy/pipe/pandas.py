from .pipe import pipe_function


@pipe_function
def df_func(df, func, *args, **kwargs):
    func_ = getattr(df, func)
    return func_(*args, **kwargs)
