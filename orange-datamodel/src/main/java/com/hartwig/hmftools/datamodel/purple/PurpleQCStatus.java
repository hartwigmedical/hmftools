package com.hartwig.hmftools.datamodel.purple;

public enum PurpleQCStatus
{
    PASS,

    WARN_DELETED_GENES,
    WARN_HIGH_COPY_NUMBER_NOISE,
    WARN_GENDER_MISMATCH,
    WARN_LOW_PURITY,

    FAIL_CONTAMINATION,
    FAIL_NO_TUMOR
}
