package com.hartwig.hmftools.esvee.common.saga;

import htsjdk.samtools.Cigar;

public record SagaMatchBySequence(
        SagaAlignment alignment
)
{
    public SagaVariant variant()
    {
        return alignment.sagaAssembly().variant();
    }

    public Cigar cigar()
    {
        return alignment.cigar();
    }
}
