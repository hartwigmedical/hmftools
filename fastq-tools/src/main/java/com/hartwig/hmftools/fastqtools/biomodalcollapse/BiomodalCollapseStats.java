package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static java.lang.String.format;

import java.util.concurrent.atomic.AtomicInteger;

public class BiomodalCollapseStats
{
    public AtomicInteger ProcessedFastqPairCount;
    public AtomicInteger WrittenConsensusFastqCount;

    public BiomodalCollapseStats()
    {
        ProcessedFastqPairCount = new AtomicInteger();
        WrittenConsensusFastqCount = new AtomicInteger();
    }

    @Override
    public String toString()
    {
        return format("BiomodalCollapse stats: processedFastqPairCount(%d) writtenConsensusFastqCount(%d)", ProcessedFastqPairCount.get(), WrittenConsensusFastqCount.get());
    }
}
