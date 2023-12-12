from sklearn.base import BaseEstimator

"""
This module contains sklearn estimators that do nothing. These are useful as 'bypass' estimators amongst a collection
of actual estimators, specifically when a polymorphic method e.g. transform() needs to be called from all estimators.
"""

class PassthroughTransformer(BaseEstimator):
    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return X

    def fit_transform(self, X, y=None):
        return self.transform(X)

    def set_output(self, transform=None):
        return self


class PassthroughClassifier(PassthroughTransformer):
    def __init__(self):
        super().__init__()

    def predict_proba(self, X, y=None):
        return X

    def predict(self, X, y=None):
        return self.predict_proba(X).max(axis=1)


class PassthroughRegression(PassthroughTransformer):
    def __init__(self):
        super().__init__()

    def predict(self, X, y=None):
        return X
