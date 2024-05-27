package com.hartwig.hmftools.redux.common;

import static java.lang.Math.round;

public class DuplicateFrequency
{
    public final int DuplicateCount;
    public long Frequency;
    public int DualStrandFrequency;

    public DuplicateFrequency(final int duplicateCount)
    {
        DuplicateCount = duplicateCount;
        Frequency = 0;
        DualStrandFrequency = 0;
    }

    public static int roundFrequency(int frequency)
    {
        if(frequency <= 10)
            return frequency;
        else if(frequency <= 100)
            return round(frequency/10) * 10;
        else if(frequency <= 1000)
            return round(frequency/100) * 100;
        else
            return round(frequency/1000) * 1000;
    }
}
