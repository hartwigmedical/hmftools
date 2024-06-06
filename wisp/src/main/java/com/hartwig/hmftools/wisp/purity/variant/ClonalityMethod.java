package com.hartwig.hmftools.wisp.purity.variant;

public enum ClonalityMethod
{
    NONE,
    VAF_PEAK,
    NO_PEAK,
    LOW_COUNT;

    public static boolean isRecomputed(final ClonalityMethod method)
    {
        switch(method)
        {
            case VAF_PEAK:
            case LOW_COUNT:
                return true;

            default:
                return false;
        }
    }

    public static boolean convertVafToPurity(final ClonalityMethod method)
    {
        return method == LOW_COUNT;
    }
}
