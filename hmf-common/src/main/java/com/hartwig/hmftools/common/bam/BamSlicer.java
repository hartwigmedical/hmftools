package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SAM_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamSlicer
{
    private final BamSlicerFilter mFilter;

    private volatile boolean mConsumerHalt = false; // allow consumer to halt processing

    public BamSlicer(int minMappingQuality)
    {
        this(minMappingQuality, false, false, false);
    }

    public BamSlicer(int minMappingQuality, boolean keepDuplicates, boolean keepSupplementaries, boolean keepSecondaries)
    {
        mFilter = new BamSlicerFilter(minMappingQuality, keepDuplicates, keepSupplementaries, keepSecondaries);
    }

    public BamSlicer(final BamSlicerFilter filter)
    {
        mFilter = filter;
    }

    public void setKeepUnmapped()
    {
        mFilter.setKeepUnmapped();
    }

    public void setKeepHardClippedSecondaries()
    {
        mFilter.setKeepHardClippedSecondaries();
    }

    public void haltProcessing() { mConsumerHalt = true; }

    public void slice(final SamReader samReader, final ChrBaseRegion region, final Consumer<SAMRecord> consumer)
    {
        slice(samReader, List.of(region), consumer);
    }

    public void slice(final SamReader samReader, final List<ChrBaseRegion> regions, final Consumer<SAMRecord> consumer)
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

                if(mFilter.passesFilters(record))
                {
                    consumer.accept(record);
                }
            }
        }
    }

    public List<SAMRecord> slice(final SamReader samReader, final ChrBaseRegion region)
    {
        return slice(samReader, createIntervals(List.of(region), samReader.getFileHeader()));
    }

    public List<SAMRecord> slice(final SamReader samReader, final QueryInterval[] queryIntervals)
    {
        final List<SAMRecord> records = Lists.newArrayList();

        if(queryIntervals == null)
            return records;

        try(final SAMRecordIterator iterator = samReader.queryOverlapping(queryIntervals))
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if(mFilter.passesFilters(record))
                {
                    records.add(record);
                }
            }
        }

        return records;
    }

    public void queryUnmapped(final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        mConsumerHalt = false;

        try(final SAMRecordIterator iterator = samReader.queryUnmapped())
        {
            while(!mConsumerHalt && iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if(mFilter.passesFilters(record))
                {
                    consumer.accept(record);
                }
            }
        }
    }

    public List<SAMRecord> queryMates(final SamReader samReader, final List<SAMRecord> records)
    {
        List<SAMRecord> mateRecords = Lists.newArrayListWithExpectedSize(records.size());

        for(SAMRecord record : records)
        {
            SAMRecord mateRecord = queryMate(samReader, record);
            if(mateRecord != null && mFilter.passesFilters(mateRecord))
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

            if(firstInPair != firstInPair(nextRecord))
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

        return mateRecord != null && mFilter.passesFilters(mateRecord) ? mateRecord : null;
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

        boolean isFirstInPair = firstInPair(record);
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

            if(isFirstInPair == firstInPair(nextRecord))
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

        return mateRecord != null && mFilter.passesFilters(mateRecord) ? mateRecord : null;
    }

    public static QueryInterval[] createIntervals(final List<ChrBaseRegion> regions, final SAMFileHeader header)
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

    // farm out the asynchronous query tasks into the ExecutorService.
    // \returns a CompletableFuture<Void> object, which the user can wait on.
    // usage:
    //    Collection<ChrBaseRegion> partitions = ...
    //
    //    BamSlicer bamSlicer = new BamSlicer(config.MinMappingQuality, false, false, false);
    //    CompletableFuture<Void> bamSliceTasks = bamSlicer.queryAsync(new File(config.BamPath), readerFactory, partitions,
    //            false, executorService, this::processRead, this::regionComplete);
    //
    //    //wait for all tasks to complete
    //    bamSliceTasks.get();
    //
    public CompletableFuture<Void> queryAsync(final File bamFile, final SamReaderFactory readerFactory, Collection<ChrBaseRegion> regions,
            boolean queryUnmapped, final ExecutorService executorService, final BiConsumer<SAMRecord, ChrBaseRegion> samRecordConsumer,
            @Nullable final Consumer<ChrBaseRegion> regionCompleteNotify)
    {
        SAM_LOGGER.info("queryAsync of {} bam regions from {}", regions.size(), bamFile.getPath());

        List<SamReader> samReaderList = Collections.synchronizedList(new ArrayList<>());

        // create one bam reader per thread
        final ThreadLocal<SamReader> threadBamReader = ThreadLocal.withInitial(() -> {
            SamReader samReader = readerFactory.open(bamFile);
            samReaderList.add(samReader);
            return samReader;
        });

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        if(queryUnmapped)
        {
            if(!mFilter.keepUnmapped())
            {
                SAM_LOGGER.warn("queryAsync on BamSlicer with keepUnmapped=false");
            }
            Runnable task = () -> {
                SAM_LOGGER.info("queryAsync unmapped");
                queryUnmapped(threadBamReader.get(), samRecord -> samRecordConsumer.accept(samRecord, null));
            };
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        for(ChrBaseRegion region : regions)
        {
            Runnable task = () -> {
                SAM_LOGGER.printf(Level.INFO, "queryAsync region(%s:%,d-%,d)", region.chromosome(), region.start(), region.end());
                slice(threadBamReader.get(), region, samRecord -> samRecordConsumer.accept(samRecord, region));
                if(regionCompleteNotify != null)
                {
                    regionCompleteNotify.accept(region);
                }
            };
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        CompletableFuture<Void> future = CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new));

        // add a task that clean up all the bam readers after all the bam reading tasks are completed
        future = future.thenRun(() -> {
            // close all sam readers
            for(SamReader samReader : samReaderList)
            {
                SAM_LOGGER.info("cleaning up");
                try
                {
                    samReader.close();
                }
                catch(IOException e)
                {
                    throw new UncheckedIOException(e);
                }
            }
            SAM_LOGGER.info("queryAsync of {} partitions complete", regions.size());
        });

        return future;
    }

    // overload for users who do not need the region complete consumer
    public CompletableFuture<Void> queryAsync(final File bamFile, final SamReaderFactory readerFactory, Collection<ChrBaseRegion> regions,
            boolean queryUnmapped, final ExecutorService executorService, final BiConsumer<SAMRecord, ChrBaseRegion> samRecordConsumer)
    {
        return queryAsync(bamFile, readerFactory, regions, queryUnmapped, executorService, samRecordConsumer, null);
    }
}
