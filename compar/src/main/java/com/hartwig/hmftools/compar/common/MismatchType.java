package com.hartwig.hmftools.compar.common;

public enum MismatchType
{
    REF_ONLY,
    NEW_ONLY,
    VALUE,
    FULL_MATCH,
    INVALID_REF, // from a missing or invalid input source
    INVALID_NEW,
    INVALID_BOTH,
    INVALID_ERROR
}
