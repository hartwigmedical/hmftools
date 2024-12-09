package com.hartwig.hmftools.redux.write;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

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
        String regionsStr = mRegions.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
        return format("%d: %d regions: %s", mRegionIndex, mRegions.size(), regionsStr);
    }
}
