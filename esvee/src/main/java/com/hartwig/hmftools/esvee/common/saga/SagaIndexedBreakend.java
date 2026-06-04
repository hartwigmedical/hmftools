package com.hartwig.hmftools.esvee.common.saga;

import com.hartwig.hmftools.common.region.BasePosition;

// For matching to SAGA variants by genomic location.
public record SagaIndexedBreakend(
        BasePosition basePosition,
        SagaVariant variant
) implements Comparable<SagaIndexedBreakend>
{
    public String chromosome()
    {
        return basePosition.Chromosome;
    }

    public int position()
    {
        return basePosition.Position;
    }

    @Override
    public int compareTo(final SagaIndexedBreakend other)
    {
        return Integer.compare(position(), other.position());
    }
}
