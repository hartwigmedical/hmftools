package com.hartwig.hmftools.cup.rna;

public class RnaCohortData
{
    public final String CancerType;
    public final int ReadLength;
    public final int SampleCount;

    public RnaCohortData(final String cancerType, final int readLength, int sampleCount)
    {
        CancerType = cancerType;
        ReadLength = readLength;
        SampleCount = sampleCount;
    }
}
