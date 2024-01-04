## Enable from <module> import *
from os.path import dirname, basename, isfile, join
import glob
modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [
    basename(f)[:-3]
    for f in modules
    if isfile(f) and not f.endswith('__init__.py')
]

## Increase pandas df display width
import pandas as pd
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 300)

## Disable interactive matplotlib
import matplotlib
matplotlib.use('Agg')

## Ignore warnings from plotting
import warnings
warnings.filterwarnings("ignore", module = "matplotlib\..*" )
warnings.filterwarnings("ignore", module = "plotnine\..*" )