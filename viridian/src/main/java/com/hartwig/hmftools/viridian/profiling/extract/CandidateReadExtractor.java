package com.hartwig.hmftools.viridian.profiling.extract;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.common.perf.TaskExecutor.executeRunnables;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.Streams.partitionStream;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_EXTRACTION_PARTITION_SIZE;

import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Consumer;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.CRAIEntry;
import htsjdk.samtools.util.FileExtensions;

// Extracts potentially viral reads from the tumor BAM and writes them to single-ended FASTA ready for alignment.
// Workers extract from disjoint regions in parallel. Each writes its own FASTA shard which is joined at the end.
public class CandidateReadExtractor
{
    @Nullable
    private final String mRefGenomeFile; // required only for CRAM decode
    private final CandidateReadFilter mFilter;
    private final int mThreads;

    private static final Logger LOGGER = LogManager.getLogger(CandidateReadExtractor.class);

    public CandidateReadExtractor(@Nullable String refGenomeFile, CandidateReadFilter filter, int threads)
    {
        mRefGenomeFile = refGenomeFile;
        mFilter = filter;
        mThreads = threads;
    }

    public void extract(String tumorBamFile, String outputFastaFile)
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
        }
        finally
        {
            workers.forEach(Worker::close);
        }
    }

    private Queue<WorkerTask> scanTasks(SamReaderFactory readerFactory, String tumorBamFile)
    {
        try(SamReader reader = readerFactory.open(new File(tumorBamFile)))
        {
            if(!reader.hasIndex())
            {
                throw new UserInputError("Tumor BAM/CRAM is not indexed: " + tumorBamFile);
            }
            // The fully unmapped block can be large, so queue it first to avoid a long tail.
            Queue<WorkerTask> tasks = new ConcurrentLinkedQueue<>(fullyUnmappedTasks(reader, tumorBamFile));
            for(SAMSequenceRecord sequence : reader.getFileHeader().getSequenceDictionary().getSequences())
            {
                if(mFilter.isCandidateContig(sequence.getSequenceName()))
                {
                    partitionChromosome(sequence, VIRAL_READ_EXTRACTION_PARTITION_SIZE)
                            .forEach(partition -> tasks.add(new RegionTask(partition)));
                }
            }
            return tasks;
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to open tumor BAM", e);
        }
    }

    private List<? extends WorkerTask> fullyUnmappedTasks(SamReader reader, String tumorBamFile)
    {
        if(reader.type() == SamReader.Type.CRAM_TYPE)
        {
            // For CRAM, the index allows efficiently querying subsections of the unmapped block, which can be read in parallel.
            return cramUnmappedPartitionTasks(tumorBamFile, mThreads);
        }
        else
        {
            // For BAM, there is no easy way to efficiently shard the unmapped block, so have to read it all sequentially.
            return List.of(new WholeUnmappedBlockTask());
        }
    }

    // Split the CRAM unmapped block into sections which can be read in parallel.
    private static List<CramUnmappedPartitionTask> cramUnmappedPartitionTasks(String cramFile, int taskCount)
    {
        File indexFile = SamFiles.findIndex(new File(cramFile));
        if(!indexFile.getName().endsWith(FileExtensions.CRAM_INDEX))
        {
            throw new UserInputError("Tumor CRAM index must be a CRAI: " + indexFile);
        }

        // Find the file offsets of the "containers" (CRAM chunks) containing unmapped reads.
        List<Long> containerOffsets;
        try(InputStream indexStream = new FileInputStream(indexFile))
        {
            containerOffsets = CRAMCRAIIndexer.readIndex(indexStream).getCRAIEntries().stream()
                    .filter(entry -> entry.getSequenceId() == NO_ALIGNMENT_REFERENCE_INDEX)
                    .map(CRAIEntry::getContainerStartByteOffset)
                    .distinct()
                    .sorted()
                    .toList();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to read tumor CRAM index: " + indexFile, e);
        }

        // Group consecutive containers into file spans. Each span then becomes a thread task.
        int containersPerSpan = max(1, containerOffsets.size() / taskCount);
        List<BAMFileSpan> spans = partitionStream(containerOffsets.stream(), containersPerSpan)
                // A span holds virtual file offsets: the container's byte offset in the top 48 bits, and an offset inside it in the low 16 bits.
                // htsjdk reads every container that starts within the span and rejects a span whose start equals its end. So the end takes
                // the last container's offset with the low 16 bits all set.
                .map(offsets -> new BAMFileSpan(new Chunk(offsets.get(0) << 16, ((offsets.get(offsets.size() - 1)) << 16) | 0xFFFF)))
                .toList();

        return IntStream.range(0, spans.size())
                .mapToObj(i -> new CramUnmappedPartitionTask(spans.get(i), i, spans.size()))
                .toList();
    }

    private interface WorkerTask
    {
        void run(SamReader reader, Consumer<SAMRecord> recordConsumer);
    }

    private record RegionTask(
            ChrBaseRegion region
    ) implements WorkerTask
    {
        @Override
        public void run(SamReader reader, Consumer<SAMRecord> recordConsumer)
        {
            try(SAMRecordIterator records = reader.queryOverlapping(region.Chromosome, region.start(), region.end()))
            {
                records.stream()
                        // If the record straddles the region boundary, then queryOverlapping will return it twice, so need to additionally
                        // filter by the alignment start position.
                        .filter(record -> record.getAlignmentStart() >= region.start())
                        .forEach(recordConsumer);
            }
        }

        @NotNull
        @Override
        public String toString() { return region.toString(); }
    }

    // For BAM. A BAM index cannot query inside the unmapped block, so the whole block is read as one task.
    private record WholeUnmappedBlockTask() implements WorkerTask
    {
        @Override
        public void run(SamReader reader, Consumer<SAMRecord> recordConsumer)
        {
            try(SAMRecordIterator records = reader.queryUnmapped())
            {
                records.forEachRemaining(recordConsumer);
            }
        }

        @NotNull
        @Override
        public String toString() { return "fully unmapped"; }
    }

    // For CRAM. One run of whole containers from the unmapped block. Mapped reads in a mixed first container are dropped,
    // because region tasks own them.
    private record CramUnmappedPartitionTask(
            BAMFileSpan containers,
            int index,
            int count
    ) implements WorkerTask
    {
        @Override
        public void run(SamReader reader, Consumer<SAMRecord> recordConsumer)
        {
            try(SAMRecordIterator records = reader.indexing().iterator(containers))
            {
                records.stream()
                        .filter(record -> record.getReferenceIndex() == NO_ALIGNMENT_REFERENCE_INDEX)
                        .forEach(recordConsumer);
            }
        }

        @NotNull
        @Override
        public String toString() { return format("fully unmapped %d/%d", index + 1, count); }
    }

    // One extraction thread's private BAM reader and FASTA shard, so scanning and writing run lock-free.
    private class Worker implements Runnable
    {
        private final TaskQueue<WorkerTask> mTasks;
        private final SamReader mReader;
        private final FastaPart mPart;
        private long mRecordCount;

        private Worker(TaskQueue<WorkerTask> tasks, SamReader reader, FastaPart part)
        {
            mTasks = tasks;
            mReader = reader;
            mPart = part;
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
            long startRecordCount = mRecordCount;
            int startCount = mPart.readCount();

            task.run(mReader, this::processRecord);

            LOGGER.debug(
                    "extract({}) {} records, {} candidate reads in {}s",
                    task, mRecordCount - startRecordCount, mPart.readCount() - startCount, format("%.1f", secondsSinceNow(startTimeMs)));
        }

        private void processRecord(SAMRecord record)
        {
            ++mRecordCount;
            if(mFilter.isCandidateRecord(record))
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

    private static int joinFastaParts(List<FastaPart> parts, String outputFastaFile)
    {
        LOGGER.debug("Joining {} candidate FASTA parts", parts.size());

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

        return candidateCount;
    }
}
