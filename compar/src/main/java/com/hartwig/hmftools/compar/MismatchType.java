package com.hartwig.hmftools.compar;

public enum MismatchType
{
    REF_ONLY,
    NEW_ONLY,
    VALUE;

    public static boolean isPresence(final MismatchType type)
    {
        return type == REF_ONLY || type == NEW_ONLY;
    }
}
