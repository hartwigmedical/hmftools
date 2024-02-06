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

    public void processRead(final SAMRecord read)
    {
        final boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();
        final boolean isSecondary = read.isSecondaryAlignment();
        final boolean isSupp = read.getSupplementaryAlignmentFlag();
        final boolean isDuplicate = read.getDuplicateReadFlag();
        final boolean isMapped = !read.getReadUnmappedFlag();
        final boolean isPaired = read.getReadPairedFlag();
        final boolean isProperPair = read.getProperPairFlag();

        increment(FlagStatType.TOTAL, passesQC);

        if(isDuplicate)
        {
            increment(FlagStatType.DUPLICATE, passesQC);
        }

        if(isMapped)
        {
            increment(FlagStatType.MAPPED, passesQC);
        }

        if(isSecondary)
        {
            increment(FlagStatType.SECONDARY, passesQC);
            return;
        }

        if(isSupp)
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

        if(!isPaired)
            return;

        // It is a paired primary read.
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

        // It is a mapped paired primary read.
        if(isProperPair)
        {
            increment(FlagStatType.PROPERLY_PAIRED, passesQC);
        }

        if(read.getMateUnmappedFlag())
        {
            increment(FlagStatType.SINGLETON, passesQC);
            return;
        }

        // It is a mapped paired primary read with a mapped mate.
        increment(FlagStatType.PAIR_MAPPED, passesQC);

        if(read.getReferenceName().equals(read.getMateReferenceName()))
            return;

        // It is a mapped paired primary read with a mapped mate toa different chromosome.
        increment(FlagStatType.INTER_CHR_PAIR_MAPPED, passesQC);

        if(read.getMappingQuality() >= MAP_QUAL_THRESHOLD)
        {
            increment(FlagStatType.INTER_CHR_PAIR_MAP_QUAL_GE5, passesQC);
        }
    }


    public void registerConsensusRead(final SAMRecord read)
    {
        final boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();
        final boolean isSupp = read.getSupplementaryAlignmentFlag();

        decrement(FlagStatType.DUPLICATE, passesQC);

        if(!isSupp)
        {
            decrement(FlagStatType.PRIMARY_DUPLICATE, passesQC);
        }
    }

    public FlagQCStats getStat(final FlagStatType type) { return mTypeStats.get(type.ordinal()); }
    public String statAsString(final FlagStatType type) { return mTypeStats.get(type.ordinal()).toString(); }
}
