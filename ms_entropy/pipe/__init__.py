from .pipe import *
from .entropy import *
from .functions import *
from .functions import print_ as print
from .pandas import df_func

# from . import df
import pandas as pd


class Object:
    def __init__(self):
        for _func in dir(pd.DataFrame):
            if _func.startswith("_"):
                continue
            setattr(self, _func, pipe_function(getattr(pd.DataFrame, _func)))

    def __getitem__(self, item):
        return pd.DataFrame.__getitem__, [item], {}


df = Object()
