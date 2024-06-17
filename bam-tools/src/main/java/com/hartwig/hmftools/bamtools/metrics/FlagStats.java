package com.hartwig.hmftools.bamtools.metrics;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMRecord;

public class FlagStats
{
    private final List<FlagQCStats> mTypeStats;

    private static final int MAP_QUAL_THRESHOLD = 5;

    public FlagStats()
    {
        mTypeStats = Arrays.stream(FlagStatType.values()).map(x -> new FlagQCStats()).collect(Collectors.toList());
    }

    public List<FlagQCStats> typeStats() { return mTypeStats; }
    public int passCount(final FlagStatType type) { return mTypeStats.get(type.ordinal()).getPassed(); }

    public void merge(final FlagStats other)
    {
        for(FlagStatType type : FlagStatType.values())
        {
            mTypeStats.get(type.ordinal()).merge(other.typeStats().get(type.ordinal()));
        }
    }

    private void increment(final FlagStatType type, boolean passes)
    {
        mTypeStats.get(type.ordinal()).record(passes);
    }

    private void decrement(final FlagStatType type, boolean passes)
    {
        mTypeStats.get(type.ordinal()).decrement(passes);
    }

    public void processRead(final SAMRecord read, boolean isConsensusRead)
    {
        boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();
        boolean isSupplementary = read.getSupplementaryAlignmentFlag();
        boolean isSecondary = read.getSupplementaryAlignmentFlag();

        if(isConsensusRead)
        {
            // consensus reads decrement the extra duplicate read as previously explained in BamReader, but
            // being extra reads they do not contribute to supplementary, mapped or other counts
            decrement(FlagStatType.DUPLICATE, passesQC);

            if(!isSupplementary && !isSecondary)
                decrement(FlagStatType.PRIMARY_DUPLICATE, passesQC);

            return;
        }

        increment(FlagStatType.TOTAL, passesQC);

        boolean isDuplicate = read.getDuplicateReadFlag();
        boolean isMapped = !read.getReadUnmappedFlag();

        if(isDuplicate)
        {
            increment(FlagStatType.DUPLICATE, passesQC); // duplicate of either primary or supplementary
        }

        if(isMapped)
        {
            increment(FlagStatType.MAPPED, passesQC);
        }

        if(read.isSecondaryAlignment())
        {
            increment(FlagStatType.SECONDARY, passesQC);
            return;
        }

        if(isSupplementary)
        {
            increment(FlagStatType.SUPPLEMENTARY, passesQC);
            return;
        }

        increment(FlagStatType.PRIMARY, passesQC);

        if(isDuplicate)
        {
            increment(FlagStatType.PRIMARY_DUPLICATE, passesQC);
        }

        if(isMapped)
        {
            increment(FlagStatType.PRIMARY_MAPPED, passesQC);
        }

        if(!read.getReadPairedFlag())
            return;

        increment(FlagStatType.PAIRED, passesQC);

        if(read.getFirstOfPairFlag())
        {
            increment(FlagStatType.READ1, passesQC);
        }

        if(read.getSecondOfPairFlag())
        {
            increment(FlagStatType.READ2, passesQC);
        }

        if(!isMapped)
            return;

        if(read.getProperPairFlag())
        {
            increment(FlagStatType.PROPERLY_PAIRED, passesQC);
        }

        if(read.getMateUnmappedFlag())
        {
            increment(FlagStatType.SINGLETON, passesQC);
            return;
        }

        increment(FlagStatType.PAIR_MAPPED, passesQC);

        if(!read.getReferenceName().equals(read.getMateReferenceName()))
        {
            increment(FlagStatType.INTER_CHR_PAIR_MAPPED, passesQC);

            if(read.getMappingQuality() >= MAP_QUAL_THRESHOLD)
            {
                increment(FlagStatType.INTER_CHR_PAIR_MAP_QUAL_GE5, passesQC);
            }
        }
    }

    public FlagQCStats getStat(final FlagStatType type) { return mTypeStats.get(type.ordinal()); }
    public String statAsString(final FlagStatType type) { return mTypeStats.get(type.ordinal()).toString(); }
}
