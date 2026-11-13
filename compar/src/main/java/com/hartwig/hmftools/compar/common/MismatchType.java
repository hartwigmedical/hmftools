package com.hartwig.hmftools.compar.common;

public enum MismatchType
{
    OLD_ONLY,
    NEW_ONLY,
    VALUE,
    FULL_MATCH,
    ABSENT_FOUND,
    ABSENT_CONFIRMED,
    INVALID_OLD, // from a missing or invalid input source
    INVALID_NEW,
    INVALID_BOTH,
    INVALID_ERROR
}
