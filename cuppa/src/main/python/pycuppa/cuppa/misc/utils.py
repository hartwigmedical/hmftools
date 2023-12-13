import pandas as pd
import numpy as np
from numpy._typing import NDArray
from typing import Iterable
import logging
logger = logging.getLogger(__name__)


def max_col(X):
    max_col_index = np.argmax(X, axis=1)
    if not isinstance(X, pd.DataFrame):
        return max_col_index
    return np.array(X.columns[max_col_index])


def get_top_cols(
    array_2d: pd.DataFrame | NDArray,
    top_n: int | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:

    ## Detect input data type
    if isinstance(array_2d, pd.DataFrame):
        index = array_2d.index
        columns = array_2d.columns
    elif isinstance(array_2d, np.ndarray):
        if len(array_2d.shape) != 2:
            logger.error("numpy array must be 2D")
            raise TypeError
        index = np.array(range(array_2d.shape[0]))
        columns = np.array(range(array_2d.shape[1]))
    else:
        logger.error("`array_2d` must be a pandas dataframe or 2D numpy array")
        raise TypeError

    ## Fill NAs
    is_na = pd.isnull(array_2d) if isinstance(array_2d, pd.DataFrame) else np.isnan(array_2d)
    array_2d = np.where(is_na, -np.inf, array_2d)

    ## Get sort order
    row_indexes = np.fliplr(np.argsort(array_2d))
    if top_n is not None:
        row_indexes = row_indexes[:, :top_n]

    ## Apply sort
    names = []
    values = []
    for i, row_indexes_i in enumerate(row_indexes):
        # i=0
        # row_indexes_i=row_indexes[i]
        names.append(columns[row_indexes_i])
        values.append(array_2d[i, row_indexes_i])

    names = np.vstack(names)
    values = np.vstack(values)
    is_minus_inf = values == -np.inf

    ## Post process
    matrices = dict(names=names, values=values)

    for key, matrix in matrices.items():
        #matrix = matrices["names"]

        ## Propagate NAs
        matrix = np.where(is_minus_inf, np.nan, matrix)

        ## Add row and column names
        matrix = pd.DataFrame(matrix, index=index)
        matrix.columns = matrix.columns + 1
        matrix.columns.name = "rank"

        matrices[key] = matrix

    names, values = matrices.values()

    return names, values


def bin_array(
    x: NDArray | pd.Series,
    bins: NDArray = np.linspace(0, 1, 11),
    right: bool = True,
    incl_outer: bool = True,
    label_precision: int = 3,
    label_type: str = "interval"
) -> pd.Categorical:
    # if False:
    #     x = [-1, 0, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 2]
    #     x = pred_summ["pred_prob_1"]
    #     bins = np.linspace(0, 1, 11)
    #     right = True
    #     incl_outer = True
    #     label_precision = 3
    #     label_type="interval"

    ## Sanitize inputs --------------------------------
    x = np.array(x)

    if len(bins)<2:
        logger.error("`bins` must be an array-like with at least 2 elements")
        raise TypeError

    if label_type not in ["interval","midpoint"]:
        logger.error("`label_type` must be 'interval' or 'midpoint'")
        raise ValueError

    bins = np.sort(np.array(bins))

    ## Make integer bins --------------------------------
    ## out of bounds values are assigned 0 when less or len(bins) when greater
    x_bins = np.array(np.digitize(x, bins, right))
    x_bins[x_bins==len(bins)] = 0 ## Also set greater out of bounds bin to 0

    if incl_outer:
        if right:
            ## Right is already included, therefore need to include left
            outer_bin_value = bins[0]
            x_bins[x == outer_bin_value] = 1
        else:
            ## Left is already included, therefore need to include right
            outer_bin_value = bins[-1]
            x_bins[x == outer_bin_value] = len(bins)-1

    ## Make bin names --------------------------------
    if label_type=="interval":
        bins_str = np.array(np.round(bins, label_precision), dtype=str)
        bin_template = "(%s,%s]" if right else "[%s,%s)"
        bin_names = [
            bin_template % (bins_str[i], bins_str[i+1])
            for i in range(len(bins)-1)
        ]

        if incl_outer:
            if right:
                bin_names[0] = "[" + bin_names[0][1:]
            else:
                bin_name = bin_names[len(bin_names)-1]
                bin_name = bin_name[:len(bin_name)-1] + ']'
                bin_names[len(bin_names)-1] = bin_name
    else:
        bin_names = [
            (bins[i] + bins[i + 1]) / 2
            for i in range(len(bins) - 1)
        ]

    bin_names = [np.nan] + bin_names
    bin_names = pd.Series(bin_names)

    ## Assign bin names --------------------------------
    x_bin_names = bin_names[x_bins]
    if label_type=="midpoint":
        return x_bin_names.values

    x_bin_names = bin_names[x_bins]
    if x_bin_names.isna().any():
        logger.warning("Could not return Categorical. `x` contains NAs or values outside of `bins`")
    else:
        x_bin_names = pd.Categorical(x_bin_names, bin_names.dropna().values)

    return x_bin_names


def check_required_columns(
    df: pd.DataFrame,
    required_columns: Iterable,
    from_index: bool = False
) -> None:

    columns_from_df = df.columns if not from_index else df.index.names

    required_columns = pd.Series(required_columns)
    missing_columns = required_columns[ ~required_columns.isin(columns_from_df) ]

    if len(missing_columns)>0:
        raise KeyError("Missing required columns: " + ", ".join(missing_columns))


def as_categorical(series: pd.Series):
    return pd.Categorical(series, categories=series.unique())