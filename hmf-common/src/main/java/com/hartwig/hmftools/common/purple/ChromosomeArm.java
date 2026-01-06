package com.hartwig.hmftools.common.purple;

public enum ChromosomeArm
{
    P_ARM,
    Q_ARM,
    CENTROMERE,
    UNKNOWN;

    public static final String CHROMOSOME_ARM_P = "P";
    public static final String CHROMOSOME_ARM_Q = "Q";

    public static String asStr(final ChromosomeArm arm)
    {
        if(arm == P_ARM)
            return CHROMOSOME_ARM_P;
        else if(arm == Q_ARM)
            return CHROMOSOME_ARM_Q;
        else if(arm == CENTROMERE)
            return CENTROMERE.toString();
        else
            return UNKNOWN.toString();
    }

    public static ChromosomeArm fromString(final String str)
    {
        if(str.equals(CHROMOSOME_ARM_P))
            return P_ARM;
        else if(str.equals(CHROMOSOME_ARM_Q))
            return Q_ARM;
        else if(str.equals(CENTROMERE.toString()))
            return CENTROMERE;
        return UNKNOWN;
    }
}
