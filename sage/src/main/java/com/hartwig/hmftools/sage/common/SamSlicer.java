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

public class SamSlicer
{
    private final List<ChrBaseRegion> mRegions;

    private final BamSlicer mBamSlicer;

    public SamSlicer(final int minMappingQuality, final ChrBaseRegion slice)
    {
        mBamSlicer = new BamSlicer(minMappingQuality, false, false, false);
        mRegions = Collections.singletonList(slice);
    }

    public SamSlicer(final int minMappingQuality, final ChrBaseRegion slice, final List<BaseRegion> panel)
    {
        mBamSlicer = new BamSlicer(minMappingQuality);
        mRegions = Lists.newArrayList();

        for(final BaseRegion panelRegion : panel)
        {
            if(panelRegion.start() <= slice.end() && panelRegion.end() >= slice.start())
            {
                ChrBaseRegion overlap = new ChrBaseRegion(slice.Chromosome,
                        Math.max(panelRegion.start(), slice.start()),
                        Math.min(panelRegion.end(), slice.end()));

                mRegions.add(overlap);
            }
        }
    }

    public void slice(final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        mBamSlicer.slice(samReader, mRegions, consumer);
    }
}
