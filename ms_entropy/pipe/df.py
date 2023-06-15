import pandas as pd
from .pipe import pipe_function


for _func in dir(pd.DataFrame):
    if _func.startswith("_"):
        continue
    globals()[_func] = pipe_function(getattr(pd.DataFrame, _func))
