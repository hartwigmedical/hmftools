package com.hartwig.hmftools.isofox.common;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.SvRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class BamSlicer
{
    private final int mMinMappingQuality;
    private boolean mKeepDuplicates;
    private boolean mKeepSupplementaries;
    private boolean mKeepSecondaries;
    private boolean mConsumerHalt; // allow consumer to halt processing

    public BamSlicer(int minMappingQuality, boolean keepDuplicates, boolean keepSupplementaries, boolean keepSecondaries)
    {
        mMinMappingQuality = minMappingQuality;
        mKeepDuplicates = keepDuplicates;
        mKeepSupplementaries = keepSupplementaries;
        mKeepSecondaries = keepSecondaries;
        mConsumerHalt = false;
    }

    public void haltProcessing() { mConsumerHalt = true; }

    public void slice(@NotNull final SamReader samReader, final List<SvRegion> regions, @NotNull final Consumer<SAMRecord> consumer)
    {
        mConsumerHalt = false;

        final QueryInterval[] queryIntervals = createIntervals(regions, samReader.getFileHeader());

        try (final SAMRecordIterator iterator = samReader.queryOverlapping(queryIntervals))
        {
            while (!mConsumerHalt && iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if (passesFilters(record))
                {
                    consumer.accept(record);
                }
            }
        }
    }

    public List<SAMRecord> slice(@NotNull final SamReader samReader, final QueryInterval[] queryIntervals)
    {
        final List<SAMRecord> records = Lists.newArrayList();

        try (final SAMRecordIterator iterator = samReader.queryOverlapping(queryIntervals))
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if (passesFilters(record))
                {
                    records.add(record);
                }
            }
        }

        return records;
    }

    private static QueryInterval[] createIntervals(@NotNull final List<SvRegion> regions, @NotNull final SAMFileHeader header)
    {
        final QueryInterval[] queryIntervals = new QueryInterval[regions.size()];

        for (int i = 0; i < regions.size(); ++i)
        {
            final SvRegion region = regions.get(i);
            int sequenceIndex = header.getSequenceIndex(region.Chromosome);

            if (sequenceIndex < 0)
                return null;


            queryIntervals[i] = new QueryInterval(sequenceIndex, region.start(), region.end());
        }

        return queryIntervals;
    }

    private boolean passesFilters(@NotNull final SAMRecord record)
    {
        if(record.getMappingQuality() < mMinMappingQuality || record.getReadUnmappedFlag())
            return false;

        if(record.isSecondaryAlignment() && !mKeepSecondaries)
            return false;

        if(record.getSupplementaryAlignmentFlag() && !mKeepSupplementaries)
            return false;

        if(record.getDuplicateReadFlag() && !mKeepDuplicates)
            return false;

        return true;
    }
}
