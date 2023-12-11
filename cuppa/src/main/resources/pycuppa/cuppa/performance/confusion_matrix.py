from __future__ import annotations

import os
from functools import cached_property
from typing import Optional

import numpy as np
import pandas as pd
import plotnine as p9

from cuppa.classifier.cuppa_prediction import CuppaPredSummary
from cuppa.constants import PAN_CANCER_CLASS_NAME
from cuppa.logger import LoggerMixin
from cuppa.misc.plotting import PlotnineFigExporter
from cuppa.performance.performance_stats import PerformanceStats


class ConfusionMatrix(LoggerMixin):
    def __init__(self, pred_summ: CuppaPredSummary, clf_name: str = None):
        """
        Confusion matrix

        A matrix where rows are predicted classes and columns are actual classes. The diagonal shows per actual class
        the correct classifications, and the off-diagonal shows the incorrect classifications.

        Parameters
        ----------
        pred_summ: CuppaPredSummary, child of pandas DataFrame
            Dataframe summarizing the top-N predictions

        clf_name: str
            Classifier name. When `pred_summ` contains the predictions from multiple classifiers, `clf_name` needs to be
            specified as a confusion matrix can only made per classifier

        Properties
        ----------
            counts_matrix: Number of correct/incorrect classifications per actual class
            props_matrix: Proportion of correct/incorrect classifications per actual class


        Methods
        -------
            plot: Plot confusion matrix as a heatmap correct/incorrect proportion as the fill colors, and the
            correct/incorrect counts as the labels.

        """

        self.clf_name = clf_name
        self.pred_summ = pred_summ

        self._select_one_clf()

    def _select_one_clf(self) -> None:
        n_clfs = self.pred_summ.clf_names.__len__()

        if n_clfs == 1:
            return None

        if self.clf_name is None and n_clfs > 1:
            self.logger.error("`pred_summ` contains results of multiple classifiers. Please specify `clf_name`")
            raise Exception

        if self.clf_name not in self.pred_summ.clf_names:
            self.logger.error("Valid values for `clf_name` are: " + ", ".join(self.pred_summ.clf_names))
            raise ValueError

        self.pred_summ = self.pred_summ[self.pred_summ["clf_name"] == self.clf_name]

    @cached_property
    def _performance(self) -> PerformanceStats:
        perf = self.pred_summ.performance()
        perf = perf[perf["clf_name"] == self.clf_name]
        return perf

    def _get_stat_all(self, column: str) -> pd.DataFrame:
        stat = self._performance.loc[self._performance["class"] == PAN_CANCER_CLASS_NAME, column]
        return stat.squeeze()

    def _make_confusion_matrix(self, normalize: bool = False) -> pd.DataFrame:
        ## Subset input data
        df = self.pred_summ[["actual_class","pred_class_1"]]

        ## Remove rows where sample does not have a prediction
        df = df[~df["pred_class_1"].isna()]

        ## Get unique classes
        classes = []
        classes += df["actual_class"].unique().tolist()
        classes += df["pred_class_1"].unique().tolist()
        classes = pd.Series(classes).sort_values().unique()

        df["actual_class"] = pd.Categorical(df["actual_class"], classes)
        df["pred_class_1"] = pd.Categorical(df["pred_class_1"], classes)

        ## Confusion matrix
        matrix = pd.crosstab(
            df["pred_class_1"], df["actual_class"],
            normalize=False if not normalize else "columns"
        )

        ## Add pan-cancer stat
        matrix.insert(0, PAN_CANCER_CLASS_NAME, np.nan) ## Add empty rows

        matrix = matrix.T ## Add empty columns
        matrix.insert(0, PAN_CANCER_CLASS_NAME, np.nan)
        matrix = matrix.T

        stat_all = self._get_stat_all("n_correct")
        if normalize:
            stat_all /= self._get_stat_all("n_total")
        matrix.iloc[0, 0] = stat_all

        if not normalize:
            matrix = matrix.astype("Int64")

        return matrix

    @cached_property
    def counts_matrix(self) -> pd.DataFrame:
        return self._make_confusion_matrix(normalize=False)

    @cached_property
    def props_matrix(self) -> pd.DataFrame:
        return self._make_confusion_matrix(normalize=True)

    def __repr__(self) -> str:
        props_and_counts = self.counts_matrix.astype(str) + " (" + self.props_matrix.round(2).astype(str) + ")"
        return props_and_counts.__repr__()

    def plot(
        self,
        path: Optional[str] = None,
        width: int = 14,
        height: int = 10,
        dpi: int = 300,
        verbose: bool = False
    ) -> None:

        if path is None:
            path = os.path.expanduser('~/Desktop/plot.pdf')

        if verbose:
            self.logger.info("Plotting confusion matrix to path: " + path)

        ## Plot data --------------------------------
        cell_values = self.props_matrix
        cell_labels = self.counts_matrix.astype(str).replace("<NA>", "")

        ## Add summary stats
        perf = self._performance.set_index("class")[["precision","recall"]].transpose()
        perf.index = [".Precision", ".Recall"]

        cell_values = pd.concat([cell_values, perf])

        cell_labels = pd.concat([
            cell_labels,
            perf.round(2).astype(str)
        ])

        ## To long form
        plot_data = pd.concat([
            cell_values.stack()._set_name("value"),
            cell_labels.stack()._set_name("label")
        ], axis=1)

        plot_data.index.names = ["row", "col"]
        plot_data.reset_index(inplace=True)

        ## Force order of rows and columns
        plot_data["row"] = pd.Categorical(plot_data["row"], np.flip(cell_values.index))
        plot_data["col"] = pd.Categorical(plot_data["col"], cell_values.columns)

        ## Plot --------------------------------
        FILL_COLORS = ['#225EA8', '#1D91C0', '#41B6C4', '#7FCDBB', '#C7E9B4', '#EDF8B1', '#FFFFD9']

        fig = (
            p9.ggplot(plot_data, p9.aes(y="row", x="col"))
            + p9.geom_tile(p9.aes(fill="value"), color="black")
            + p9.geom_text(p9.aes(label="label"), size=7)
            + p9.scale_fill_gradientn(colors=FILL_COLORS, limits=[0, 1])
            + p9.scale_x_discrete(expand=(0, 0.5), name="Actual class")
            + p9.scale_y_discrete(expand=(0, 0.5), name="Predicted class")
            + p9.theme_bw()
            + p9.theme(
                panel_grid=p9.element_blank(),
                axis_text_x=p9.element_text(angle=270, vjust=1),
                legend_position="none"
            )
        ).draw(show=False, return_ggplot=False)

        ## Export
        fig_exporter = PlotnineFigExporter(width=width, height=height, dpi=dpi)
        fig_exporter.export(fig, path)

