package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.BiConsumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AsyncBamLociReader
{
    private static final Logger logger = LogManager.getLogger(AsyncBamLociReader.class);

    // use inheritance to save some memory, it is a bit dirty
    // but Amber is running out of memory
    static class Task<E extends GenomePosition> extends ArrayList<E>
    {
        String chromosome;
        int positionStart;
        int positionEnd;
    }

    static class BamReaderThread<E extends GenomePosition> extends Thread
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

                try (final SAMRecordIterator iterator = mSamReader.queryOverlapping(task.chromosome, task.positionStart,
                        task.positionEnd))
                {
                    while (iterator.hasNext())
                    {
                        final SAMRecord record = iterator.next();

                        if (!passesFilters(record))
                        {
                            continue;
                        }

                        int alignmentEnd = record.getAlignmentEnd();

                        for (E genomePosition : task)
                        {
                            // check overlap
                            if (genomePosition.position() >= record.getAlignmentStart() && genomePosition.position() <= alignmentEnd)
                            {
                                mConsumer.accept(genomePosition, record);
                            }
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

            if(record.isSecondaryAlignment())
                return false;

            if(record.getSupplementaryAlignmentFlag())
                return false;

            return !record.getDuplicateReadFlag();
        }
    }

    public static <E extends GenomePosition> void processBam(final String bamFile, final SamReaderFactory samReaderFactory,
            final List<E> loci, BiConsumer<E, SAMRecord> asyncRecordHandler, int threadCount, int minMappingQuality)
            throws InterruptedException
    {
        logger.debug("Processing {} potential sites in bam {}", loci.size(), bamFile);

        final Queue<Task<E>> taskQ = new ConcurrentLinkedQueue<>();
        
        // create genome regions from the loci
        populateTaskQueue(loci, taskQ);

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

        AmberTaskCompletion taskCompletion = new AmberTaskCompletion(taskQ.size());
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
    public static <E extends GenomePosition> void populateTaskQueue(final List<E> sortedGenomePositions, final Queue<Task<E>> taskQ)
    {
        int minGap = 2000;
        List<Task<E>> tasks = new ArrayList<>();

        while (true)
        {
            tasks.clear();
            Task<E> task = new Task<>();
            task.chromosome = null;
            task.positionStart = -1;
            task.positionEnd = -1;

            for (E position : sortedGenomePositions)
            {
                if (task.chromosome == null || !task.chromosome.equals(position.chromosome()) || task.positionEnd + minGap < position.position())
                {
                    if (task.chromosome != null)
                    {
                        // we finalise previous one
                        tasks.add(task);
                        task = new Task<>();
                    }

                    task.chromosome = position.chromosome();
                    task.positionStart = position.position();
                }

                if (task.positionEnd > position.position())
                {
                    // this means the input genome positions are not sorted
                    throw new RuntimeException("Genome position going backwards, input might not be sorted");
                }

                // add to existing region
                task.positionEnd = position.position();
                task.add(position);
            }

            if (task.chromosome != null)
            {
                // we finalise previous one
                tasks.add(task);
            }

            if (tasks.size() <= AmberConstants.OPTIMAL_BAM_SLICE_REGIONS)
            {
                break;
            }
            minGap += 200;
        }
        taskQ.addAll(tasks);

        logger.info("{} loci, {} genome regions, min gap = {}", sortedGenomePositions.size(), taskQ.size(), minGap);
    }
}
