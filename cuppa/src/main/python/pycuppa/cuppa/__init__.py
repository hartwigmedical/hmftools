## Enable from <module> import *
from os.path import dirname, basename, isfile, join
import glob
modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [
    basename(f)[:-3]
    for f in modules
    if isfile(f) and not f.endswith('__init__.py')
]

## Set up logging
from cuppa.logger import initialize_logging
initialize_logging()

## Increase pandas df display width
import pandas as pd
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 300)