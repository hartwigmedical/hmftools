package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.common.perf.TaskExecutor.executeRunnables;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.virusdetect.VirusConstants.EXTRACTION_PARTITION_SIZE;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Extracts potentially viral reads from the tumor BAM and writes them to single-ended FASTA ready for alignment.
// The scan is split into jobs which workers claim as they finish, each worker writing its own FASTA shard.
public class CandidateReadExtractor
{
    @Nullable
    private final String mRefGenomeFile; // required only for CRAM decode
    private final CandidateReadFilter mFilter;
    private final int mThreads;

    private static final Logger LOGGER = LogManager.getLogger(CandidateReadExtractor.class);

    public CandidateReadExtractor(@Nullable String refGenomeFile, CandidateReadFilter filter)
    {
        this(refGenomeFile, filter, 1);
    }

    public CandidateReadExtractor(@Nullable String refGenomeFile, CandidateReadFilter filter, int threads)
    {
        mRefGenomeFile = refGenomeFile;
        mFilter = filter;
        mThreads = threads;
    }

    public int extractToFasta(String tumorBamFile, String outputFastaFile)
    {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        if(mRefGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSequence(new File(mRefGenomeFile));
        }

        List<Consumer<Worker>> jobs = scanJobs(readerFactory, tumorBamFile);
        List<Worker> workers = createWorkers(min(mThreads, jobs.size()), readerFactory, tumorBamFile, outputFastaFile);

        int candidateCount;
        try
        {
            runJobs(jobs, workers);
            candidateCount = joinFastaParts(workers, outputFastaFile);
        }
        finally
        {
            workers.forEach(Worker::close);
        }

        LOGGER.info("Extracted {} candidate reads to {}", candidateCount, outputFastaFile);
        return candidateCount;
    }

    private List<Consumer<Worker>> scanJobs(SamReaderFactory readerFactory, String tumorBamFile)
    {
        try(SamReader reader = readerFactory.open(new File(tumorBamFile)))
        {
            if(reader.hasIndex())
            {
                return shardedScanJobs(reader.getFileHeader().getSequenceDictionary());
            }

            if(mThreads > 1)
            {
                throw new UserInputError("Multi-threaded extraction requires an indexed BAM/CRAM: " + tumorBamFile);
            }
            return List.of(Worker::scanAll);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to open tumor BAM", e);
        }
    }

    private static List<Consumer<Worker>> shardedScanJobs(SAMSequenceDictionary dictionary)
    {
        List<Consumer<Worker>> jobs = new ArrayList<>();

        // The unmapped reads sit in one block scanned by a single long-running job. Claimed first so it runs alongside
        // the region jobs from the start, rather than tacking its full duration onto the end.
        jobs.add(Worker::scanUnmapped);

        for(SAMSequenceRecord sequence : dictionary.getSequences())
        {
            for(ChrBaseRegion region : partitionChromosome(sequence, EXTRACTION_PARTITION_SIZE))
            {
                jobs.add(worker -> worker.scanRegion(region));
            }
        }
        return jobs;
    }

    private List<Worker> createWorkers(int workerCount, SamReaderFactory readerFactory, String tumorBamFile, String outputFastaFile)
    {
        BamSlicer slicer = new BamSlicer(0, true, true, true);
        slicer.setKeepUnmapped();
        Predicate<SAMRecord> isCandidate = record ->
                !record.getDuplicateReadFlag() && !record.isSecondaryOrSupplementary() && mFilter.isCandidate(record);

        List<Worker> workers = new ArrayList<>();
        try
        {
            for(int i = 0; i < workerCount; ++i)
            {
                workers.add(new Worker(
                        readerFactory.open(new File(tumorBamFile)), FastaPart.create(outputFastaFile, i), slicer, isCandidate));
            }
        }
        catch(RuntimeException e)
        {
            workers.forEach(Worker::close);
            throw e;
        }
        return workers;
    }

    private static void runJobs(List<Consumer<Worker>> jobs, List<Worker> workers)
    {
        AtomicInteger nextJob = new AtomicInteger();
        List<Runnable> tasks = workers.stream()
                .map(worker -> (Runnable) () ->
                {
                    int index;
                    while((index = nextJob.getAndIncrement()) < jobs.size())
                    {
                        jobs.get(index).accept(worker);
                    }
                })
                .toList();

        if(!executeRunnables(tasks, workers.size()))
        {
            throw new RuntimeException("Candidate read extraction failed");
        }
    }

    private static int joinFastaParts(List<Worker> workers, String outputFastaFile)
    {
        long startTimeMs = System.currentTimeMillis();
        int candidateCount = 0;
        try(OutputStream out = new BufferedOutputStream(new FileOutputStream(outputFastaFile)))
        {
            for(Worker worker : workers)
            {
                FastaPart part = worker.part();
                part.close();
                Files.copy(part.path(), out);
                candidateCount += part.readCount();
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to join candidate FASTA parts", e);
        }

        LOGGER.debug("Joined {} FASTA parts in {}s", workers.size(), format("%.1f", secondsSinceNow(startTimeMs)));
        return candidateCount;
    }

    // One extraction thread's private BAM reader and FASTA shard, so scanning and writing run lock-free.
    private static class Worker
    {
        private final SamReader mReader;
        private final FastaPart mPart;
        private final BamSlicer mSlicer;
        private final Predicate<SAMRecord> mIsCandidate;

        private Worker(SamReader reader, FastaPart part, BamSlicer slicer, Predicate<SAMRecord> isCandidate)
        {
            mReader = reader;
            mPart = part;
            mSlicer = slicer;
            mIsCandidate = isCandidate;
        }

        private FastaPart part()
        {
            return mPart;
        }

        // Whole-file linear scan, the only option without an index.
        private void scanAll()
        {
            for(SAMRecord record : mReader)
            {
                addIfCandidate(record);
            }
        }

        private void scanRegion(ChrBaseRegion region)
        {
            long startTimeMs = System.currentTimeMillis();
            int startCount = mPart.readCount();

            // A mapped read is owned by the partition containing its start, so copies returned by an overlapping
            // neighbour partition are ignored.
            mSlicer.slice(mReader, region, record ->
            {
                if(record.getAlignmentStart() >= region.start())
                {
                    addIfCandidate(record);
                }
            });

            LOGGER.debug("region({}) {} candidates in {}s",
                    region, mPart.readCount() - startCount, format("%.1f", secondsSinceNow(startTimeMs)));
        }

        private void scanUnmapped()
        {
            long startTimeMs = System.currentTimeMillis();
            int startCount = mPart.readCount();

            mSlicer.queryUnmapped(mReader, this::addIfCandidate);

            LOGGER.debug("Unmapped reads {} candidates in {}s",
                    mPart.readCount() - startCount, format("%.1f", secondsSinceNow(startTimeMs)));
        }

        private void addIfCandidate(SAMRecord record)
        {
            if(mIsCandidate.test(record))
            {
                mPart.add(record);
            }
        }

        private void close()
        {
            mPart.discard();
            try
            {
                mReader.close();
            }
            catch(IOException e)
            {
                LOGGER.warn("failed to close tumor BAM: {}", e.getMessage());
            }
        }
    }
}
