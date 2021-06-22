package com.hartwig.hmftools.sage.sam;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SamSlicer
{
    private final List<BaseRegion> mRegions;

    private final BamSlicer mBamSlicer;

    public SamSlicer(final int minMappingQuality, final BaseRegion slice)
    {
        mBamSlicer = new BamSlicer(minMappingQuality);
        mRegions = Collections.singletonList(slice);
    }

    public SamSlicer(final int minMappingQuality, final BaseRegion slice, final List<BaseRegion> panel)
    {
        mBamSlicer = new BamSlicer(minMappingQuality);
        mRegions = Lists.newArrayList();

        for(final BaseRegion panelRegion : panel)
        {
            if(slice.Chromosome.equals(panelRegion.Chromosome) && panelRegion.start() <= slice.end()
                    && panelRegion.end() >= slice.start())
            {
                BaseRegion overlap = new BaseRegion(slice.Chromosome,
                        Math.max(panelRegion.start(), slice.start()),
                        Math.min(panelRegion.end(), slice.end()));

                mRegions.add(overlap);
            }
        }
    }

    public void slice(final SamReader samReader, @NotNull final Consumer<SAMRecord> consumer)
    {
        mBamSlicer.slice(samReader, mRegions, consumer);

        // TODO: is this required or are regions already exclusive??
        // return QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));


        /*
        final QueryInterval[] queryIntervals = createIntervals(mRegions, samReader.getFileHeader());

        try(final SAMRecordIterator iterator = samReader.queryOverlapping(queryIntervals))
        {
            while(iterator.hasNext())
            {
                final SAMRecord record = iterator.next();
                if(samRecordMeetsQualityRequirements(record))
                {
                    consumer.accept(record);
                }
            }
        }
        */
    }
}
