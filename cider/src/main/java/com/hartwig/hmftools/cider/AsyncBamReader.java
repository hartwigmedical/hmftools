package com.hartwig.hmftools.cider;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.BiConsumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AsyncBamReader
{
    private static final Logger logger = LogManager.getLogger(AsyncBamReader.class);

    // use inheritance to save some memory, it is a bit dirty
    // but Amber is running out of memory
    static class Task<E extends GenomeRegion>
    {
        E genomeRegion;
    }

    static class BamReaderThread<E extends GenomeRegion> extends Thread
    {
        final Queue<Task<E>> mTaskQ;
        final SamReader mSamReader;
        final int mMinMappingQuality;
        BiConsumer<E, SAMRecord> mConsumer;

        BamReaderThread(final String bamFile, final SamReaderFactory samReaderFactory, final Queue<Task<E>> inTaskQ,
                int minMappingQuality, BiConsumer<E, SAMRecord> consumer)
        {
            mTaskQ = inTaskQ;
            mSamReader = samReaderFactory.open(new File(bamFile));
            mMinMappingQuality = minMappingQuality;
            mConsumer = consumer;
        }

        @Override
        public void run()
        {
            logger.debug("bam reader thread start");

            while (true)
            {
                Task<E> task;
                try
                {
                    task = mTaskQ.remove();
                }
                catch (NoSuchElementException e)
                {
                    // finished processing
                    break;
                }

                try (final SAMRecordIterator iterator = task.genomeRegion == null ?
                        mSamReader.queryUnmapped() :
                        mSamReader.queryOverlapping(task.genomeRegion.chromosome(), task.genomeRegion.start(), task.genomeRegion.end()))
                {
                    while (iterator.hasNext())
                    {
                        final SAMRecord record = iterator.next();

                        if (!passesFilters(record))
                        {
                            continue;
                        }

                        //if (task.genomeRegion.start() >= record.getAlignmentStart() && task.genomeRegion.end() <= alignmentEnd)
                        if (task.genomeRegion.start() <= record.getAlignmentEnd() && task.genomeRegion.end() >= record.getAlignmentStart())
                        {
                            mConsumer.accept(task.genomeRegion, record);
                        }
                    }
                }
            }

            try {
                mSamReader.close();
            }
            catch (IOException e)
            {
                logger.error("IO exception in SamReader::close: {}", e.getMessage());
            }

            logger.debug("bam reader thread finish");
        }

        private boolean passesFilters(final SAMRecord record)
        {
            if(record.getMappingQuality() < mMinMappingQuality || record.getReadUnmappedFlag())
                return false;

            return !record.getDuplicateReadFlag();
        }
    }

    public static <E extends GenomeRegion> void processBam(final String bamFile, final SamReaderFactory samReaderFactory,
            final Collection<E> regions, BiConsumer<E, SAMRecord> asyncRecordHandler, int threadCount, int minMappingQuality)
            throws InterruptedException
    {
        logger.debug("Processing {} potential sites in bam {}", regions.size(), bamFile);

        final Queue<Task<E>> taskQ = new ConcurrentLinkedQueue<>();
        
        // create genome regions from the loci
        populateTaskQueue(regions, taskQ);
        BamTaskCompletion taskCompletion = new BamTaskCompletion(taskQ.size());

        // we create the consumer and producer
        var bamReaders = new ArrayList<BamReaderThread<E>>();

        for (int i = 0; i < Math.max(threadCount, 1); ++i)
        {
            var t = new BamReaderThread<>(bamFile, samReaderFactory, taskQ, minMappingQuality, asyncRecordHandler);
            t.setName(String.format("worker-%d", i));
            t.start();
            bamReaders.add(t);
        }

        logger.info("{} bam reader threads started", bamReaders.size());

        for (BamReaderThread<E> t : bamReaders)
        {
            while (t.isAlive())
            {
                // check status every 30 seconds
                t.join(30_000);

                // check status
                taskCompletion.progress(taskQ.size());
            }
        }

        logger.info("{} bam reader threads finished", bamReaders.size());
    }

    /**
     * Group the genome regions
     */
    public static <E extends GenomeRegion> void populateTaskQueue(final Collection<E> genomeRegions, final Queue<Task<E>> taskQ)
    {
        for (E genomeRegion : genomeRegions)
        {
            Task<E> task = new Task<>();
            task.genomeRegion = genomeRegion;
            taskQ.add(task);
        }
        // another one for unmapped reads
        taskQ.add(new Task<>());
    }
}
