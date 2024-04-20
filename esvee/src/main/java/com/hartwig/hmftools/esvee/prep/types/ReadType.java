package com.hartwig.hmftools.esvee.prep.types;

public enum ReadType
{
    NO_SUPPORT,
    CANDIDATE_SUPPORT,
    SUPPORT,
    EXACT_SUPPORT,
    JUNCTION,
    EXPECTED,
    RECOVERED;

    public static int rank(final ReadType type)
    {
        switch(type)
        {
            case JUNCTION: return 7;
            case EXACT_SUPPORT: return 6;
            case SUPPORT: return 5;
            case CANDIDATE_SUPPORT: return 4;
            case NO_SUPPORT: return 0;
            default: return 0;
        }

    }
}
