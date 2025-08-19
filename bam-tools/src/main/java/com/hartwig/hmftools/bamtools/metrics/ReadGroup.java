package com.hartwig.hmftools.bamtools.metrics;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class ReadGroup
{
    public final boolean IsConsensus;
    public final int MaxReadStart;
    public final List<int[]> CombinedAlignedBaseCoords;

    public ReadGroup(final int maxReadStart, boolean isConsensus)
    {
        IsConsensus = isConsensus;
        CombinedAlignedBaseCoords = Lists.newArrayListWithCapacity(2);
        MaxReadStart = maxReadStart;
    }
}
