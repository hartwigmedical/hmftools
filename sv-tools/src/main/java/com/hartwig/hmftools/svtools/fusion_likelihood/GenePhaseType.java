package com.hartwig.hmftools.svtools.fusion_likelihood;

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

    public static GenePhaseType intAsType(int phase)
    {
        switch(phase)
        {
            case 0: return PHASE_NON_CODING;
            case 1 : return PHASE_5P_UTR;
            case 2: return PHASE_0;
            case 3: return PHASE_1;
            case 4: return PHASE_2;
        }

        return PHASE_5P_UTR;
    }

    public static boolean isPhasedType(GenePhaseType type)
    {
        return type == PHASE_0 || type == PHASE_1 || type == PHASE_2;
    }
}
