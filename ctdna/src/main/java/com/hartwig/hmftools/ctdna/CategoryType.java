package com.hartwig.hmftools.ctdna;

public enum CategoryType
{
    REFERENCE,
    FUSION,
    GERMLINE_SV,
    AMP_DEL,
    GERMLINE_MUTATION,
    REPORTABLE_MUTATION,
    OTHER_SV,
    SUBCLONAL_MUTATION,
    OTHER_CODING_MUTATION,
    OTHER_MUTATION;

    public boolean isSV() { return this == FUSION || this == AMP_DEL || this == OTHER_SV; }

    public boolean isMutation() { return this == REPORTABLE_MUTATION || this == OTHER_MUTATION || this == OTHER_CODING_MUTATION
            || this == SUBCLONAL_MUTATION; }
}
