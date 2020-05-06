package com.hartwig.hmftools.isofox.common;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class BamSlicer
{
    private final int mMinMappingQuality;
    private boolean mDropDuplicates;
    private boolean mDropSupplementaries;
    private boolean mConsumerHalt; // allow consumer to halt processing

    public BamSlicer(int minMappingQuality, boolean dropDuplicates, boolean dropSecondSupps)
    {
        mMinMappingQuality = minMappingQuality;
        mDropDuplicates = dropDuplicates;
        mDropSupplementaries = dropSecondSupps;
        mConsumerHalt = false;
    }

    public void haltProcessing() { mConsumerHalt = true; }

    public void slice(@NotNull final SamReader samReader, final List<GenomeRegion> regions, @NotNull final Consumer<SAMRecord> consumer)
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

    public void slice(@NotNull final SamReader samReader, final QueryInterval[] queryIntervals, @NotNull final Consumer<SAMRecord> consumer)
    {
        mConsumerHalt = false;

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

    private static QueryInterval[] createIntervals(@NotNull final List<GenomeRegion> regions, @NotNull final SAMFileHeader header)
    {
        final QueryInterval[] queryIntervals = new QueryInterval[regions.size()];

        for (int i = 0; i < regions.size(); ++i)
        {
            final GenomeRegion region = regions.get(i);
            int sequenceIndex = header.getSequenceIndex(region.chromosome());

            if (sequenceIndex < 0)
                return null;


            queryIntervals[i] = new QueryInterval(sequenceIndex, (int) region.start(), (int) region.end());
        }

        return queryIntervals;
    }

    private boolean passesFilters(@NotNull final SAMRecord record)
    {
        if(record.getMappingQuality() < mMinMappingQuality || record.getReadUnmappedFlag())
            return false;

        if(mDropSupplementaries && (record.getSupplementaryAlignmentFlag() || record.isSecondaryAlignment()))
            return false;

        if(mDropDuplicates && record.getDuplicateReadFlag())
            return false;

        return true;
    }
}
