package com.hartwig.hmftools.bachelor.types;

public enum FilterType
{
    NONE,
    PASS,
    VCF_FILTERED,
    ARTEFACT,
    BLACKLIST,
    CLINVAR_CONFLICTING,
    CLINVAR_BENIGN,
    CLINVAR_LIKELY_BENIGN,
    CLINVAR_UNANNOTATED;
}
