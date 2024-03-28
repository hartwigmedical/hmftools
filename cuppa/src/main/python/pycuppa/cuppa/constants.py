from __future__ import annotations

try:
    from importlib import resources as impresources
except ImportError:
    ## Try `importlib_resources` backported for Python<3.7
    import importlib_resources as impresources

## Classifier ================================
class SUB_CLF_NAMES:
    GEN_POS = "gen_pos"
    SNV96 = "snv96"
    EVENT = "event"

    GENE_EXP = "gene_exp"
    ALT_SJ = "alt_sj"

    @classmethod
    def get_dna_clf_names(cls) -> list[str]:
        return [cls.GEN_POS, cls.SNV96, cls.EVENT]

    @classmethod
    def get_rna_clf_names(cls) -> list[str]:
        return [cls.GENE_EXP, cls.ALT_SJ]

SIG_QUANTILE_TRANSFORMER_NAME = "sig"

class META_CLF_NAMES:
    DNA_COMBINED = "dna_combined"
    RNA_COMBINED = "rna_combined"
    COMBINED = "combined"


class LAYER_NAMES:
    SUB_CLF = "sub_clfs"
    META_CLF = "meta_clfs"
    COMBINED = "prob_combiner"


class CLF_GROUPS:
    COMBINED = "combined"
    DNA = "dna"
    RNA = "rna"

    @classmethod
    def get_all(cls) -> list[str]:
        return [cls.COMBINED, cls.DNA, cls.RNA]

    MAPPINGS_CLF_NAMES = {
        META_CLF_NAMES.COMBINED: COMBINED,
        META_CLF_NAMES.DNA_COMBINED: DNA,
        META_CLF_NAMES.RNA_COMBINED: RNA,

        SUB_CLF_NAMES.GEN_POS: DNA,
        SUB_CLF_NAMES.SNV96: DNA,
        SUB_CLF_NAMES.EVENT: DNA,

        SUB_CLF_NAMES.GENE_EXP: RNA,
        SUB_CLF_NAMES.ALT_SJ: RNA
    }


CUPPA_PREDICTION_INDEX_NAMES = ["sample_id", "data_type", "clf_group", "clf_name", "feat_name", "feat_value"]


## Paths ================================
RESOURCES_DIR = impresources.files("cuppa")/"resources"

DEFAULT_FUSION_OVERRIDES_PATH = str(RESOURCES_DIR/"feature_overrides/20230731-fusion_overrides.split_adenoid_salivary.txt")
MOCK_DATA_DIR = str(RESOURCES_DIR/"mock_data")

RSCRIPT_PLOT_PREDICTIONS_PATH = str(impresources.files("cuppa")/"visualization/plot_predictions.R")
RSCRIPT_PLOT_CONFUSION_PATH = str(impresources.files("cuppa")/"performance/plot_confusion.R")


## Misc ================================
DEFAULT_FEATURE_PREFIX_SEPERATOR = "__"

DEFAULT_KEYWORDS_MALE_CLASSES = ("prostate","penis","penile","testis","testicular")
DEFAULT_KEYWORDS_FEMALE_CLASSES = ("ovary","fallopian","ovarian","uterus","uterine","endometrium","endometrial","cervix","cervical")
SEX_FEATURE_NAME = "event.trait.is_male"

PAN_CANCER_CLASS_NAME = ".All"

NA_FILL_VALUE = -0.00000001

