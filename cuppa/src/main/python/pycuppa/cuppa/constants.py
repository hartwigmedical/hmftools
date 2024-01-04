from __future__ import annotations

import importlib

## Classifier ================================
class SUB_CLF_NAMES:
    GEN_POS = "gen_pos"
    SNV96 = "snv96"
    EVENT = "event"

    GENE_EXP = "gene_exp"
    ALT_SJ = "alt_sj"


class META_CLF_NAMES:
    DNA_COMBINED = "dna_combined"
    RNA_COMBINED = "rna_combined"
    COMBINED = "combined"


class LAYER_NAMES:
    SUB_CLF = "sub_clfs"
    META_CLF = "meta_clfs"
    COMBINED = "prob_combiner"


CLF_GROUPS = {
    META_CLF_NAMES.COMBINED: "combined",
    META_CLF_NAMES.DNA_COMBINED: "dna",
    META_CLF_NAMES.RNA_COMBINED: "rna",

    SUB_CLF_NAMES.GEN_POS: "dna",
    SUB_CLF_NAMES.SNV96: "dna",
    SUB_CLF_NAMES.EVENT: "dna",

    SUB_CLF_NAMES.GENE_EXP: "rna",
    SUB_CLF_NAMES.ALT_SJ: "rna"
}

CLF_NAME_ALIASES = {
    SUB_CLF_NAMES.GEN_POS: "rmd"
}

CUPPA_PREDICTION_INDEX_NAMES = ["sample_id", "data_type", "clf_group", "clf_name", "feat_name", "feat_value"]


## Paths ================================
RESOURCES_DIR = importlib.resources.files("resources")

DEFAULT_CUPPA_CLASSIFIER_PATH = str(RESOURCES_DIR/"cuppa_classifier.pickle.gz")
DEFAULT_FUSION_OVERRIDES_PATH = str(RESOURCES_DIR/"feature_overrides/20230731-fusion_overrides.split_adenoid_salivary.txt")
MOCK_DATA_DIR = str(RESOURCES_DIR/"mock_data")
RSCRIPT_PLOT_PREDICTIONS_PATH = str(importlib.resources.files("cuppa")/"visualization/plot_predictions.R")


## Misc ================================
DEFAULT_FEATURE_PREFIX_SEPERATOR = "__"

DEFAULT_KEYWORDS_MALE_CLASSES = ("prostate","penis","penile","testis","testicular")
DEFAULT_KEYWORDS_FEMALE_CLASSES = ("ovary","fallopian","ovarian","uterus","uterine","endometrium","endometrial","cervix","cervical")
SEX_FEATURE_NAME = "event.trait.is_male"

PAN_CANCER_CLASS_NAME = ".All"

