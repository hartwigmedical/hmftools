package com.hartwig.hmftools.common.segmentation;

public enum Arm
{
    P,
    Q;

    public static Arm fromString(String s)
    {
        if(s.equals(P.name()))
        {
            return P;
        }
        else if(s.equals(Q.name()))
        {
            return Q;
        }
        throw new IllegalArgumentException("Unknown Arm: " + s);
    }
}
