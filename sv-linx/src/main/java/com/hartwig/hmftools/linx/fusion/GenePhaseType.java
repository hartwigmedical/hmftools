package com.hartwig.hmftools.linx.fusion;

public enum GenePhaseType
{
    PHASE_NON_CODING,
    PHASE_5P_UTR,
    PHASE_0,
    PHASE_1,
    PHASE_2;

    public static final int PHASE_MAX = 5;

    public static int typeAsInt(GenePhaseType type)
    {
        switch(type)
        {
            case PHASE_NON_CODING: return 0;
            case PHASE_5P_UTR: return 1;
            case PHASE_0: return 2;
            case PHASE_1: return 3;
            case PHASE_2: return 4;
        }

        return 0;
    }
}
