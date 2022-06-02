package com.hartwig.hmftools.sage.common;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class SamSlicer implements SamSlicerInterface
{
    private final SamReader mSamReader;
    private final List<ChrBaseRegion> mRegions;

    private final BamSlicer mBamSlicer;

    public SamSlicer(final SamReader samReader, final int minMappingQuality, final List<ChrBaseRegion> regions)
    {
        mSamReader = samReader;
        mBamSlicer = new BamSlicer(minMappingQuality);
        mRegions = regions;
    }

    /*
    public void slice(final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        mBamSlicer.slice(samReader, mRegions, consumer);
    }
    */

    public void slice(final Consumer<SAMRecord> consumer)
    {
        mBamSlicer.slice(mSamReader, mRegions, consumer);
    }
}
