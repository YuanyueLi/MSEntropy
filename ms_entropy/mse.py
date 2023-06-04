class Pipe:
    def __init__(self, df=None):
        self.dataframe = df

    def __rrshift__(self, other):
        return Pipe(other)

    def __rshift__(self, other):
        func, args, kwargs = other
        return func(self.dataframe, *args, **kwargs)


def function_decorator_dataframe(func):
    def wrapper(*args, **kwargs):
        return func, args, kwargs

    return wrapper


if __name__ == "__main__":

    @function_decorator_dataframe
    def f_test(a, b=None, c=None):
        print(a, b, c)

    "abc" >> Pipe() >> f_test("b", c="c")
