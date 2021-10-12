package com.hartwig.hmftools.pave;

public enum SpliceImpactType
{
    UNKNOWN,
    HOMOLOGY_SHIFT,
    BASE_SHIFT,
    BASE_CHANGE,
    REGION_DELETED,
    OUTSIDE_RANGE;

    public boolean isDisruptive()
    {
        switch(this)
        {
            case BASE_CHANGE:
            case BASE_SHIFT:
            case REGION_DELETED:
                return true;

            case HOMOLOGY_SHIFT:
                return false;

            case OUTSIDE_RANGE:
            default:
                return false;

        }
    }
}
