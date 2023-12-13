from __future__ import annotations

from typing import Optional, TYPE_CHECKING

import pandas as pd
import plotnine as p9

from cuppa.constants import META_CLF_NAMES, SUB_CLF_NAMES
from cuppa.logger import LoggerMixin
from cuppa.misc.plotting import PlotnineFigExporter
from cuppa.misc.utils import as_categorical

if TYPE_CHECKING:
    from cuppa.classifier.cuppa_classifier import CuppaClassifier


class FeatureImportance(pd.DataFrame, LoggerMixin):

    def __init__(self, importances: pd.DataFrame):
        super().__init__(importances)

    @property
    def _constructor(self):
        return FeatureImportance

    @staticmethod
    def from_cuppa_classifier(cuppa_classifier: CuppaClassifier, na_fill: float = 0.0):

        clfs = {
            SUB_CLF_NAMES.GEN_POS: cuppa_classifier.gen_pos_clf,
            SUB_CLF_NAMES.SNV96: cuppa_classifier.snv96_clf,
            SUB_CLF_NAMES.EVENT: cuppa_classifier.event_clf,
            SUB_CLF_NAMES.GENE_EXP: cuppa_classifier.gene_exp_clf,
            SUB_CLF_NAMES.ALT_SJ: cuppa_classifier.alt_sj_clf,

            META_CLF_NAMES.DNA_COMBINED: cuppa_classifier.dna_combined_clf,
            META_CLF_NAMES.RNA_COMBINED: cuppa_classifier.rna_combined_clf,
        }

        importances = {
            clf_name: clf["logistic_regression"].feat_imp(mode="coef")
            for clf_name, clf in clfs.items()
        }
        importances = pd.concat(importances, axis=1)

        ## NAs occur because some RNA classifiers may be missing cancer types due to lack of samples
        importances.fillna(na_fill, inplace=True)

        importances.index.name = "class"
        importances.columns.names = ["clf_name","feat_name"]

        return FeatureImportance(importances)

    def to_tsv(self, path: str) -> None:
        self.to_csv(path, sep='\t')

    @staticmethod
    def from_tsv(path: str):
        #path = "/Users/lnguyen/Desktop/feat_imp.tsv.gz"
        df = pd.read_table(path, header=[0,1], index_col=0)
        return df

    @property
    def clf_names(self):
        return self.columns.get_level_values("clf_name").unique()

    def summarize(
        self,
        clf_names: Optional[str|list[str]]= None,
        top_n: int = 15
    ) -> pd.DataFrame:

        ## Subset by clf --------------------------------
        df = self
        if clf_names is not None:
            clf_names = pd.Series(clf_names)

            if ~clf_names.isin(self.clf_names).all():
                self.logger.error("Valid values for `clf_name` are: " + ", ".join(self.clf_names))

            df = self[clf_names]

        df = df.groupby("class").agg(["mean", "std"])

        ## Reformat df --------------------------------
        ## The above groupby creates an unnamed column multi-index level
        column_index_names = list(df.columns.names)
        column_index_names[-1] = "metric"
        df.columns.names = column_index_names

        ## Long format
        df = df.reset_index().melt(id_vars="class")

        ## Force ordering
        df["clf_name"] = as_categorical(df["clf_name"])
        df["feat_name"] = as_categorical(df["feat_name"])

        ## Metric values as separate columns
        index_cols = list(df.columns[~df.columns.isin(["metric","value"])])
        df = df.pivot(index=index_cols, columns="metric")
        df.columns = df.columns.get_level_values("metric")
        df = df.reset_index()

        ## Select top n features --------------------------------
        df["mean_abs"] = df["mean"].abs()
        df["rank"] = df.groupby("class")["mean_abs"].rank(ascending=False, method="first")
        #df = df.sort_values(["class", "rank"])
        df = df[df["rank"] <= top_n]

        return df

    @staticmethod
    def get_feat_affixes(feat_names: pd.Series, sep = '[.]|__', strip_event_prefix: bool = True) -> (pd.Series, pd.Series):

        if strip_event_prefix:
            ## Use second prefix for event classifier since the feature names are like this:
            ## event.virus.HPV
            ## event.sv.LINE
            feat_names = feat_names.str.replace("event.", "", regex=False)

        affixes = feat_names.str.split(sep, regex=True, n=1)

        prefixes = [list_i[0] for list_i in affixes]
        suffixes = [list_i[1] for list_i in affixes]

        return prefixes, suffixes

    def plot(
        self,
        path: Optional[str] = None,
        clf_names: Optional[str] = None,

        width: int | float = 20,
        height: int | float = 15,
        dpi: int = 300,

        use_abs_values: bool = True,
        facet_ncol: int = 6,
        label_va: str = "bottom",

    ):
        summ = self.summarize(clf_names=clf_names)

        ## Plot data --------------------------------
        ## Get feat name affixes
        prefixes, suffixes = self.get_feat_affixes(
            summ["feat_name"],
            strip_event_prefix = clf_names is not None
        )

        summ["feat_prefix"] = prefixes
        summ["feat_suffix"] = suffixes

        ## Fills and labels
        summ["label"] = summ["feat_suffix"]
        if not use_abs_values:
            summ["value"] = summ["mean"]
        else:
            summ["value"] = summ["mean_abs"]

            summ["sign_string"] = ""
            summ.loc[summ["mean"] < 0, "sign_string"] = "(-) "
            summ["label"] = summ["sign_string"] + summ["label"]

        ## Label y position
        y_range = summ["value"].max() - summ["value"].min()
        label_ypos = summ["value"].min() + y_range * 0.01

        ##
        #summ["clf_name"]
        summ["feat_prefix"] = pd.Categorical(summ["feat_prefix"], summ["feat_prefix"].unique())

        ## Plot --------------------------------
        fig = (
            p9.ggplot(summ, p9.aes(x="rank", y="value"))
            + p9.facet_wrap("class", ncol=facet_ncol)

            + p9.geom_bar(p9.aes(fill="feat_prefix"), stat="identity", width=1, color="black", size=0.2)
            + p9.geom_linerange(p9.aes(ymin="value - std", ymax="value + std"), color="grey")
            + p9.geom_text(p9.aes(label="label", y=label_ypos), angle=90, va=label_va, size=6)

            + p9.scale_fill_brewer(type="qual", palette="Set3")
            + p9.labs(x="Feature rank", y="Coefficient (mean, SD)", fill="Feature group")
            + p9.theme_bw()
            + p9.theme(
                panel_grid=p9.element_blank(),
            )
        ).draw(show=False, return_ggplot=False)

        ## Export
        fig_exporter = PlotnineFigExporter(width=width, height=height, dpi=dpi)
        fig_exporter.export(fig, path)
