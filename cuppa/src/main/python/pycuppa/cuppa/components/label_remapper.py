import numpy as np
import pandas as pd
import warnings

from sklearn.base import BaseEstimator

## Remap labels to merge super types ================================
# if False:
#     mappings = [
#         ("Breast", "Breast triple negative"),
#         ("Lung", "Lung: Non-small Cell"),
#         ("Lung", "Lung: Small Cell"),
#         ("NET", "Lung: NET"),
#         ("NET", "Pancreas: NET"),
#         ("NET", "Small intestine/Colorectum: NET"),
#     ]
#
#     probs = pd.read_table("/Users/lnguyen/Hartwig/pycuppa/data/models/Hartwig_PCAWG/23-all_dna/09-feat_contrib/01-baseline/tables/probs.txt", index_col=0)
#     label_remapper = LabelRemapper(mappings=mappings)
#     X_trans, y_trans = label_remapper.transform(probs, y)

def _check_mappings(mappings):
    if not isinstance(mappings, list) and not isinstance(mappings[0], tuple):
        raise TypeError("`mappings` must be a list of tuples")

    if not all([len(tuple_i)==2 for tuple_i in mappings]):
        raise TypeError("Each tuple in `mappings` must be of length 2")

def remap_y_labels(y, mappings):
    _check_mappings(mappings)

    y_trans = np.array(y) ## Create a numpy array copy to prevent the original data being modified
    for new_class, old_class in mappings:
        is_old_class = y_trans==old_class
        if not any(is_old_class):
            warnings.warn("'%s' is not an existing class" % old_class)

        y_trans[is_old_class] = new_class

    if isinstance(y, pd.Series):
        y_trans = pd.Series(y_trans, y.index)

    return y_trans

def merge_X_columns(X, mappings, copy=True):
    if False:
        X = probs
        copy=True
        mappings = label_mappings

    if not isinstance(X, pd.DataFrame):
        raise ValueError("`X` must be a pandas dataframe")

    X_trans = X.copy() if copy else X
    X_trans.columns = remap_y_labels(X_trans.columns, mappings)
    X_trans = X_trans.groupby(level=0, axis=1).sum(min_count=1) ## Min count = min number of valid values. Do not do row sum if all values are NAs

    return X_trans

class LabelRemapper(BaseEstimator):
    def __init__(self, mappings, bypass=False):
        self.mappings = mappings
        self.bypass = bypass

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        X_trans = X if self.bypass else merge_X_columns(X, self.mappings)
        if y is None:
            return X_trans

        y_trans = y if self.bypass else remap_y_labels(y, self.mappings)

        return X_trans, y_trans

    def fit_transform(self, X, y=None):
        return self.transform(X, y)

    def set_output(self, transform=None):
        return self
