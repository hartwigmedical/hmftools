package com.hartwig.hmftools.esvee.common.saga;

import htsjdk.samtools.Cigar;

public record SagaMatchBySequence(
        SagaVariant variant,
        Cigar cigar,
        int alignScore
)
{
}
