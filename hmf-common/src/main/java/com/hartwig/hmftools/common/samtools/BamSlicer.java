package com.hartwig.hmftools.common.samtools;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SAM_LOGGER;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class BamSlicer
{
    private final int mMinMappingQuality;
    private final boolean mKeepDuplicates;
    private final boolean mKeepSupplementaries;
    private final boolean mKeepSecondaries;
    private boolean mKeepUnmapped;

    private volatile boolean mConsumerHalt = false; // allow consumer to halt processing

    public BamSlicer(int minMappingQuality)
    {
        this(minMappingQuality, false, false, false);
    }

    public BamSlicer(int minMappingQuality, boolean keepDuplicates, boolean keepSupplementaries, boolean keepSecondaries)
    {
        mMinMappingQuality = minMappingQuality;
        mKeepDuplicates = keepDuplicates;
        mKeepSupplementaries = keepSupplementaries;
        mKeepSecondaries = keepSecondaries;
        mKeepUnmapped = false;
    }

    public void setKeepUnmapped() { mKeepUnmapped = true; }

    public void haltProcessing() { mConsumerHalt = true; }

    public void slice(@NotNull final SamReader samReader, final List<ChrBaseRegion> regions, @NotNull final Consumer<SAMRecord> consumer)
    {
        mConsumerHalt = false;

        final QueryInterval[] queryIntervals = createIntervals(regions, samReader.getFileHeader());

        if(queryIntervals == null)
            return;

        try(final SAMRecordIterator iterator = samReader.queryOverlapping(queryIntervals))
        {
            while(!mConsumerHalt && iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if(passesFilters(record))
                {
                    consumer.accept(record);
                }
            }
        }
    }

    public List<SAMRecord> slice(@NotNull final SamReader samReader, final ChrBaseRegion region)
    {
        return slice(samReader, createIntervals(Lists.newArrayList(region), samReader.getFileHeader()));
    }

    public List<SAMRecord> slice(@NotNull final SamReader samReader, final GenomePosition variantRegion)
    {
        int position = variantRegion.position();

        final QueryInterval[] queryIntervals = createIntervals(Lists.newArrayList(
                new ChrBaseRegion(variantRegion.chromosome(), position, position)), samReader.getFileHeader());

        return slice(samReader, queryIntervals);
    }

    public List<SAMRecord> slice(@NotNull final SamReader samReader, final QueryInterval[] queryIntervals)
    {
        final List<SAMRecord> records = Lists.newArrayList();

        if(queryIntervals == null)
            return records;

        try(final SAMRecordIterator iterator = samReader.queryOverlapping(queryIntervals))
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if(passesFilters(record))
                {
                    records.add(record);
                }
            }
        }

        return records;
    }

    public List<SAMRecord> queryMates(final SamReader samReader, final List<SAMRecord> records)
    {
        List<SAMRecord> mateRecords = Lists.newArrayListWithExpectedSize(records.size());

        for(SAMRecord record : records)
        {
            SAMRecord mateRecord = queryMate(samReader, record);
            if(mateRecord != null && passesFilters(mateRecord))
                mateRecords.add(mateRecord);
        }

        return mateRecords;
    }

    public SAMRecord findRead(
            final SamReader samReader, final String readId, final String chromosome, int alignmentStart,
            boolean firstInPair, boolean supplementary, int maxReadDepth)
    {
        SAMRecordIterator iter = samReader.queryAlignmentStart(chromosome, alignmentStart);

        int readCount = 0;
        SAMRecord mateRecord = null;
        while(iter.hasNext())
        {
            SAMRecord nextRecord = iter.next();

            ++readCount;

            if(maxReadDepth > 0 && readCount >= maxReadDepth)
                break;

            if(firstInPair != nextRecord.getFirstOfPairFlag())
                continue;

            // must match supplementary status so as not to be confused with the mate of its supplementary pair
            if(nextRecord.getSupplementaryAlignmentFlag() != supplementary)
                continue;

            if(nextRecord.getReadName().equals(readId))
            {
                mateRecord = nextRecord;
                break;
            }
        }

        iter.close();

        return mateRecord != null && passesFilters(mateRecord) ? mateRecord : null;
    }

    public SAMRecord queryMate(final SamReader samReader, final SAMRecord record)
    {
        // the SAM-tools implementation is nothing special and can crash if it encounters multiple reads with the same name (see issue #1164)
        // so just implement this manually

        SAMRecordIterator iter;
        if(record.getMateReferenceIndex() == -1)
        {
            iter = samReader.queryUnmapped();
        }
        else
        {
            iter = samReader.queryAlignmentStart(record.getMateReferenceName(), record.getMateAlignmentStart());
        }

        boolean isFirstInPair = record.getFirstOfPairFlag();
        boolean isFirstSupplementary = record.getSupplementaryAlignmentFlag();

        SAMRecord mateRecord = null;
        while(iter.hasNext())
        {
            SAMRecord nextRecord = iter.next();
            if(!nextRecord.getReadPairedFlag())
            {
                if(record.getReadName().equals(nextRecord.getReadName()))
                {
                    SAM_LOGGER.error("read({}) loc({}:{}) isFirst({}) mate not paired",
                            record.getReadName(), record.getContig(), record.getAlignmentStart(), isFirstInPair);
                    return null;
                }

                continue;
            }

            if(isFirstInPair == nextRecord.getFirstOfPairFlag())
                continue;

            // must match supplementary status so as not to be confused with the mate of its supplementary pair
            if(nextRecord.getSupplementaryAlignmentFlag() != isFirstSupplementary)
                continue;

            if(record.getReadName().equals(nextRecord.getReadName()))
            {
                mateRecord = nextRecord;
                break;
            }
        }

        iter.close();

        return mateRecord != null && passesFilters(mateRecord) ? mateRecord : null;
    }

    private static QueryInterval[] createIntervals(final List<ChrBaseRegion> regions, final SAMFileHeader header)
    {
        final QueryInterval[] queryIntervals = new QueryInterval[regions.size()];

        for(int i = 0; i < regions.size(); ++i)
        {
            final ChrBaseRegion region = regions.get(i);
            int sequenceIndex = header.getSequenceIndex(region.Chromosome);

            if(sequenceIndex < 0)
            {
                SAM_LOGGER.error("cannot find sequence index for chromosome {} in bam header", region.Chromosome);
                return null;
            }

            queryIntervals[i] = new QueryInterval(sequenceIndex, region.start(), region.end());
        }

        return queryIntervals;
    }

    private boolean passesFilters(final SAMRecord record)
    {
        if(record.getMappingQuality() < mMinMappingQuality)
            return false;

        if(record.getReadUnmappedFlag() && !mKeepUnmapped)
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
