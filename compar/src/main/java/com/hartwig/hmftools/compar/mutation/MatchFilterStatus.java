package com.hartwig.hmftools.compar.mutation;

public enum MatchFilterStatus
{
    BOTH_UNFILTERED(true),
    REF_FILTERED(false),
    NEW_FILTERED(false);

    private final boolean mCanComparePurpleFields;

    MatchFilterStatus(final boolean canComparePurpleFields)
    {
        mCanComparePurpleFields = canComparePurpleFields;
    }

    public boolean canComparePurpleFields()
    {
        return mCanComparePurpleFields;
    }
}
