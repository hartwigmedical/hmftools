package com.hartwig.hmftools.redux.write;

import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PartitionInfo
{
    private final int mRegionIndex;
    private final List<ChrBaseRegion> mRegions;
    private final BamWriter mBamWriter;

    public PartitionInfo(final int index, final List<ChrBaseRegion> regions, final BamWriter bamWriter)
    {
        mRegionIndex = index;
        mRegions = regions;
        mBamWriter = bamWriter;
    }

    public int regionIndex() { return mRegionIndex; }
    public List<ChrBaseRegion> regions() { return mRegions; }
    public BamWriter bamWriter() { return mBamWriter; }

    public String toString()
    {
        return format("%d: %d regions: %s", mRegionIndex, mRegions.size(), partitionInfoStr(mRegions));
    }

    public static String partitionInfoStr(final List<ChrBaseRegion> regions)
    {
        if(regions.size() <= 10)
            return regions.stream().map(x -> x.toString()).collect(Collectors.joining(";"));

        return format("%s;%s;multiple", regions.get(0), regions.get(1));
    }

    public static boolean isAltRegionContig(final String regionContig)
    {
        return !HumanChromosome.contains(regionContig) && !MitochondrialChromosome.contains(regionContig);
    }
}
