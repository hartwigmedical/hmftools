package com.hartwig.hmftools.viridian.detection.read_extract;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.common.perf.TaskExecutor.executeRunnables;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_EXTRACTION_PARTITION_SIZE;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Extracts potentially viral reads from the tumor BAM and writes them to single-ended FASTA ready for alignment.
// Workers extract from disjoint regions in parallel. Each writes its own FASTA shard which is joined at the end.
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

        Queue<WorkerTask> tasks = scanTasks(readerFactory, tumorBamFile);
        TaskQueue<WorkerTask> queue = new TaskQueue<>(tasks, "scan tasks", 0);
        List<Worker> workers = new ArrayList<>();

        try
        {
            for(int i = 0; i < min(mThreads, tasks.size()); ++i)
            {
                workers.add(new Worker(queue, readerFactory.open(new File(tumorBamFile)), FastaPart.create(outputFastaFile, i)));
            }

            boolean execSuccess = executeRunnables(workers, workers.size());
            if(!execSuccess)
            {
                throw new RuntimeException("Candidate read extraction failed");
            }

            int candidateCount = joinFastaParts(workers.stream().map(Worker::part).toList(), outputFastaFile);
            LOGGER.debug("Extracted {} candidate reads", candidateCount);
            return candidateCount;
        }
        finally
        {
            workers.forEach(Worker::close);
        }
    }

    private static Queue<WorkerTask> scanTasks(SamReaderFactory readerFactory, String tumorBamFile)
    {
        try(SamReader reader = readerFactory.open(new File(tumorBamFile)))
        {
            if(!reader.hasIndex())
            {
                throw new UserInputError("Tumor BAM/CRAM is not indexed: " + tumorBamFile);
            }

            Queue<WorkerTask> tasks = new ConcurrentLinkedQueue<>();
            // The unmapped block is not sharded, and it can be quite large, so run it first to avoid a long tail.
            tasks.add(WorkerTask.UNMAPPED_READS);
            for(SAMSequenceRecord sequence : reader.getFileHeader().getSequenceDictionary().getSequences())
            {
                partitionChromosome(sequence, VIRAL_READ_EXTRACTION_PARTITION_SIZE)
                        .forEach(partition -> tasks.add(new WorkerTask(partition)));
            }
            return tasks;
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to open tumor BAM", e);
        }
    }

    private static int joinFastaParts(List<FastaPart> parts, String outputFastaFile)
    {
        long startTimeMs = System.currentTimeMillis();
        int candidateCount = 0;
        try(OutputStream out = new BufferedOutputStream(new FileOutputStream(outputFastaFile)))
        {
            for(FastaPart part : parts)
            {
                part.close();
                Files.copy(part.path(), out);
                candidateCount += part.readCount();
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to join candidate FASTA parts", e);
        }

        LOGGER.debug("Joined {} FASTA parts in {}s", parts.size(), format("%.1f", secondsSinceNow(startTimeMs)));
        return candidateCount;
    }

    private record WorkerTask(
            @Nullable ChrBaseRegion region
    )
    {
        static final WorkerTask UNMAPPED_READS = new WorkerTask(null);

        @NotNull
        @Override
        public String toString()
        {
            return region != null ? region.toString() : "unmapped";
        }
    }

    // One extraction thread's private BAM reader and FASTA shard, so scanning and writing run lock-free.
    private class Worker implements Runnable
    {
        private final TaskQueue<WorkerTask> mTasks;
        private final SamReader mReader;
        private final FastaPart mPart;
        private final BamSlicer mSlicer;

        private Worker(TaskQueue<WorkerTask> tasks, SamReader reader, FastaPart part)
        {
            mTasks = tasks;
            mReader = reader;
            mPart = part;
            // Ignore duplicates, supplementaries, and secondaries. These would be dropped by our filter anyway, so may
            // as well filter upfront.
            mSlicer = new BamSlicer(0, false, false, false);
            mSlicer.setKeepUnmapped();
        }

        private FastaPart part() { return mPart; }

        @Override
        public void run()
        {
            while(true)
            {
                try
                {
                    processTask(mTasks.removeItem());
                }
                catch(NoSuchElementException e)
                {
                    break;
                }
            }
        }

        private void processTask(WorkerTask task)
        {
            long startTimeMs = System.currentTimeMillis();
            int startCount = mPart.readCount();

            ChrBaseRegion region = task.region();
            if(region == null)
            {
                mSlicer.queryUnmapped(mReader, this::processRecord);
            }
            else
            {
                // A mapped read is owned by the partition containing its start, so copies returned by an overlapping
                // neighbour partition are ignored.
                mSlicer.slice(
                        mReader, region, record ->
                        {
                            if(record.getAlignmentStart() >= region.start())
                            {
                                processRecord(record);
                            }
                        });
            }

            LOGGER.debug(
                    "scan({}) {} candidates in {}s",
                    task, mPart.readCount() - startCount, format("%.1f", secondsSinceNow(startTimeMs)));
        }

        private void processRecord(SAMRecord record)
        {
            if(mFilter.isCandidate(record))
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
                LOGGER.warn("Failed to close tumor BAM: {}", e.getMessage());
            }
        }
    }
}
