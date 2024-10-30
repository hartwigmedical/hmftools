package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.common.sv.SvUtils.isDiscordant;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class PartitionStats
{
    public long TotalReads;
    public long DuplicateReads;
    public long ChimericReads;
    public long InterPartition;
    public long UnmappedReads;
    public double ProcessTime;

    public PartitionStats()
    {
        TotalReads = 0;
        DuplicateReads = 0;
        ChimericReads = 0;
        InterPartition = 0;
        UnmappedReads = 0;
        ProcessTime = 0;
    }

    public void processRead(final SAMRecord read, final ChrBaseRegion partition, boolean isConsensusRead)
    {
        if(isConsensusRead)
        {
            --DuplicateReads;
        }
        else
        {
            ++TotalReads;

            if(read.getDuplicateReadFlag())
                ++DuplicateReads;
        }

        if(read.getReadUnmappedFlag())
        {
            ++UnmappedReads;
        }
        else if(isDiscordant(read))
        {
            ++ChimericReads;

            if(!partition.containsPosition(read.getMateReferenceName(), read.getMateAlignmentStart()))
                ++InterPartition;
        }
    }
}
