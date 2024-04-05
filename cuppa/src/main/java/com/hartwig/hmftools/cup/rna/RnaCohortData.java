package com.hartwig.hmftools.cup.rna;

@Deprecated
public class RnaCohortData
{
    public final String CancerType;
    public final int SampleCount;

    public RnaCohortData(final String cancerType, int sampleCount)
    {
        CancerType = cancerType;
        SampleCount = sampleCount;
    }
}
