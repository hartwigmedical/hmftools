package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static java.lang.String.format;

import java.util.concurrent.atomic.AtomicInteger;

public class BiomodalCollapseStats
{
    public AtomicInteger ProcessedFastqPairCount;
    public AtomicInteger WrittenConsensusFastqCount;

    // TODO: remove
    public AtomicInteger BadWrittenFastqPairCount;

    public BiomodalCollapseStats()
    {
        ProcessedFastqPairCount = new AtomicInteger();
        WrittenConsensusFastqCount = new AtomicInteger();
        BadWrittenFastqPairCount = new AtomicInteger();
    }

    @Override
    public String toString()
    {
        return format("BiomodalCollapse stats: processedFastqPairCount(%d) writtenConsensusFastqCount(%d) badWrittenFastqPairCount(%d)", ProcessedFastqPairCount.get(), WrittenConsensusFastqCount.get(), BadWrittenFastqPairCount.get());
    }
}
