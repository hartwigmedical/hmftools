import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator

class RandomResampler(BaseEstimator):
    def __init__(self, up_target=None, up_ratio_thres=None, down_target=None, down_ratio_thres=None, seed=0):
        self.up_target = up_target
        self.up_ratio_thres = up_ratio_thres

        self.down_target = down_target
        self.down_ratio_thres = down_ratio_thres

        self.seed = seed

    def _calc_resample_targets(self, y, up_target=None, up_ratio_thres=None, down_target=None, down_ratio_thres=None):
        if False:
            up_target = None
            up_ratio_thres = None
            down_target = 200
            down_ratio_thres = 0.5

        classes, counts = np.unique(y, return_counts=True)
        df = pd.DataFrame({
            "class": classes,
            "count": counts,
            "direction": 0,
            "up_target": 0,
            "up_ratio": 0,
            "down_target": 0,
            "down_ratio": 0,
        })

        if up_target is not None:
            df.loc[df["count"] < up_target, "direction"] = 1

            df["up_target"] = np.where(df["direction"] == 1, up_target, 0)
            df["up_ratio"] = df["up_target"] / df["count"]

            if up_ratio_thres is not None:
                df["up_ratio"] = np.minimum(df["up_ratio"], up_ratio_thres)

            df["up_target"] = np.round(df["up_ratio"] * df["count"])
            df["up_target"] = df["up_target"].astype(int)

        if down_target is not None:
            df.loc[df["count"] > down_target, "direction"] = -1

            df["down_target"] = np.where(df["direction"] == -1, down_target, 0)
            df["down_ratio"] = df["down_target"] / df["count"]

            if down_ratio_thres is not None:
                df["down_ratio"] = np.maximum(df["down_ratio"], down_ratio_thres)

            df.loc[df["direction"] != -1, "down_ratio"] = 0
            df["down_target"] = np.round(df["down_ratio"] * df["count"])
            df["down_target"] = df["down_target"].astype(int)

        return df

    def _resample_classes(self, resample_info, y, seed=0):
        # if False:
        #     seed=0
        #     resample_info = _calc_resample_targets(y, up_target=100, up_ratio_thres=5, down_target=500, down_ratio_thres=0.5)

        if seed is not None:
            np.random.seed(seed)

        ##
        resample_targets = resample_info["up_target"] + resample_info["down_target"]
        resample_targets.index = resample_info["class"]

        ##
        row_info = pd.DataFrame({
            # "index": range(len(y)),
            "class": y,
            "name": y.index
        })
        row_info = row_info.reset_index(drop=True)

        ## Get row indexes
        row_info_new = []
        for class_i in np.unique(y):
            # class_i="CNS: Pilocytic astrocytoma"
            samples = row_info[row_info["class"] == class_i]
            # row_info_i = row_info_i.drop("class", axis=1)

            resample_target = resample_targets[class_i]

            if resample_target == 0:
                row_info_new.append(samples)
                continue

            replace = True if resample_target > samples.shape[0] else False
            indexes = np.random.choice(samples.index, size=resample_target, replace=replace)

            samples_new = samples.loc[indexes]  # .sort_index()

            samples_new["dup_suffix"] = samples_new.groupby("name").cumcount().astype(str)
            samples_new["dup_suffix"] = "_" + samples_new["dup_suffix"]
            samples_new["dup_suffix"] = samples_new["dup_suffix"].replace("_0", "")

            samples_new["name"] = samples_new["name"] + samples_new["dup_suffix"]
            samples_new = samples_new.drop("dup_suffix", axis=1)

            row_info_new.append(samples_new)

        row_info_new = pd.concat(row_info_new)
        # row_info_new = row_info_new.sort_index()

        ## Select rows
        y_trans = y[row_info_new.index]
        y_trans.index = row_info_new["name"]
        return y_trans

        # if X is not None:
        #     X_trans = X.iloc[row_info_new.index]
        #     X_trans.index = row_info_new["name"]
        #     return X_trans, y_trans
        #
        # return y_trans

    def fit_resample(self, y, X=None) -> pd.Series | tuple[pd.DataFrame, pd.Series]:
        resample_info = self._calc_resample_targets(
            y,
            up_target=self.up_target,
            up_ratio_thres=self.up_ratio_thres,
            down_target=self.down_target,
            down_ratio_thres=self.down_ratio_thres
        )
        self.resample_info = resample_info

        y_trans = self._resample_classes(resample_info=resample_info, y=y, seed=self.seed)
        X_trans = X.loc[y_trans.index]

        return X_trans, y_trans

    def fit(self, X=None, y=None):
        return self

    def transform(self, X, y=None):
        return X

    def set_output(self, transform=None):
        return self
