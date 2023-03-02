package com.hartwig.hmftools.datamodel.variant;

public interface AllelicDepth
{
    int totalReadCount();

    int alleleReadCount();

    default double alleleFrequency()
    {
        return (double) alleleReadCount() / totalReadCount();
    }
}
