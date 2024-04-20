package com.hartwig.hmftools.sage.common;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class SamSlicer implements SamSlicerInterface
{
    private final SamReader mSamReader;
    private final List<ChrBaseRegion> mRegions;

    private final BamSlicer mBamSlicer;

    public SamSlicer(final SamReader samReader, final int minMappingQuality, final List<ChrBaseRegion> regions, boolean keepSupplementaries)
    {
        mSamReader = samReader;
        mBamSlicer = new BamSlicer(minMappingQuality, false, keepSupplementaries, false);
        mRegions = regions;
    }

    public void slice(final Consumer<SAMRecord> consumer)
    {
        mBamSlicer.slice(mSamReader, mRegions, consumer);
    }
}
