package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.virusdetect.VirusConstants.EXTRACTION_PARTITION_SIZE;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

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
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        if(mRefGenomeFile != null)
        {
            factory = factory.referenceSequence(new File(mRefGenomeFile));
        }

        int candidateCount = mThreads > 1
                ? extractParallel(factory, tumorBamFile, outputFastaFile)
                : extractSequential(factory, tumorBamFile, outputFastaFile);

        LOGGER.info("Extracted {} candidate reads to {}", candidateCount, outputFastaFile);
        return candidateCount;
    }

    private int extractSequential(SamReaderFactory factory, String tumorBamFile, String outputFastaFile)
    {
        int candidateCount = 0;
        try(SamReader reader = factory.open(new File(tumorBamFile));
                BufferedWriter writer = createBufferedWriter(outputFastaFile))
        {
            for(SAMRecord record : reader)
            {
                if(isCandidate(record))
                {
                    writeFastaRecord(writer, record);
                    ++candidateCount;
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to extract candidate reads", e);
        }

        return candidateCount;
    }

    private int extractParallel(SamReaderFactory readerFactory, String tumorBamFile, String outputFastaFile)
    {
        // TODO: separate sequence dict read into a function. could even include the partition generation too
        SAMSequenceDictionary dictionary;
        try(SamReader reader = readerFactory.open(new File(tumorBamFile)))
        {
            if(!reader.hasIndex())
            {
                throw new UserInputError("Multi-threaded extraction requires an indexed BAM/CRAM: " + tumorBamFile);
            }
            dictionary = reader.getFileHeader().getSequenceDictionary();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to open tumor BAM", e);
        }

        List<ChrBaseRegion> partitions = new ArrayList<>();
        for(SAMSequenceRecord sequence : dictionary.getSequences())
        {
            partitions.addAll(partitionChromosome(sequence, EXTRACTION_PARTITION_SIZE));
        }

        BamSlicer slicer = new BamSlicer(0, true, true, true);
        slicer.setKeepUnmapped();

        // One reader and one FASTA part per worker thread, both created lazily on first use, so scanning and writing
        // run lock-free.
        List<SamReader> readers = Collections.synchronizedList(new ArrayList<>());
        List<FastaPart> parts = Collections.synchronizedList(new ArrayList<>());
        AtomicInteger nextPartIndex = new AtomicInteger();

        ThreadLocal<SamReader> threadSamReader = ThreadLocal.withInitial(() ->
        {
            SamReader reader = readerFactory.open(new File(tumorBamFile));
            readers.add(reader);
            return reader;
        });
        ThreadLocal<FastaPart> threadFastaPart = ThreadLocal.withInitial(() ->
        {
            FastaPart part = FastaPart.create(outputFastaFile, nextPartIndex.getAndIncrement());
            parts.add(part);
            return part;
        });

        ExecutorService executor = Executors.newFixedThreadPool(mThreads);
        List<CompletableFuture<Void>> futures = new ArrayList<>();

        // The unmapped reads sit in one block scanned by a single long-running task. Submitted first so it
        // runs alongside the region tasks from the start, rather than tacking its full duration onto the end.
        futures.add(CompletableFuture.runAsync(() -> readAndProcessUnmapped(slicer, threadSamReader, threadFastaPart), executor));
        for(ChrBaseRegion region : partitions)
        {
            futures.add(CompletableFuture.runAsync(() -> readAndProcessRegion(slicer, threadSamReader, threadFastaPart, region), executor));
        }

        try
        {
            CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
        }
        catch(InterruptedException | ExecutionException e)
        {
            throw new RuntimeException("Candidate extraction failed", e);
        }
        finally
        {
            executor.shutdown();
        }

        // FIXME: needs to be closed also if the above try/catch fails?
        closeReaders(readers);
        return concatFastaParts(parts, outputFastaFile);
    }

    private void readAndProcessRegion(BamSlicer slicer, ThreadLocal<SamReader> threadReader, ThreadLocal<FastaPart> threadPart,
            ChrBaseRegion region)
    {
        long startTimeMs = System.currentTimeMillis();
        FastaPart fastaPart = threadPart.get();
        int startCount = fastaPart.readCount();

        slicer.slice(threadReader.get(), region, record -> processRegionRecord(record, region, fastaPart));

        LOGGER.debug("region({}) {} candidates in {}s", region, fastaPart.readCount() - startCount, secondsSince(startTimeMs));
    }

    private void readAndProcessUnmapped(BamSlicer slicer, ThreadLocal<SamReader> threadReader, ThreadLocal<FastaPart> threadPart)
    {
        long startTimeMs = System.currentTimeMillis();
        FastaPart part = threadPart.get();
        int startCount = part.readCount();

        slicer.queryUnmapped(
                threadReader.get(), record ->
                {
                    if(isCandidate(record))
                    {
                        part.add(record);
                    }
                });

        LOGGER.debug("Unmapped reads {} candidates in {}s", part.readCount() - startCount, secondsSince(startTimeMs));
    }

    private void processRegionRecord(SAMRecord record, ChrBaseRegion region, FastaPart part)
    {
        // A mapped read is owned by the partition containing its start, so copies returned by an overlapping neighbour
        // partition are ignored.
        if(record.getAlignmentStart() >= region.start() && isCandidate(record))
        {
            part.add(record);
        }
    }

    private boolean isCandidate(SAMRecord record)
    {
        return !isExcluded(record) && mFilter.isCandidate(record);
    }

    private static boolean isExcluded(SAMRecord record)
    {
        return record.getDuplicateReadFlag() || record.isSecondaryOrSupplementary();
    }

    private static void writeFastaRecord(BufferedWriter writer, SAMRecord record) throws IOException
    {
        writer.write(">" + makeFastaLabel(record));
        writer.newLine();
        writer.write(record.getReadString());
        writer.newLine();
    }

    private static String makeFastaLabel(SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            return record.getReadName() + (record.getFirstOfPairFlag() ? "/1" : "/2");
        }
        else
        {
            return record.getReadName();
        }
    }

    // TODO: put in general utils class
    private static String secondsSince(long startTimeMs)
    {
        return String.format("%.1f", (System.currentTimeMillis() - startTimeMs) / 1000.0);
    }

    // TODO: worth putting in a separate file for readability?
    // A worker's private FASTA shard, written lock-free by its owning thread and concatenated byte-wise into the
    // output. Not thread-safe.
    private static class FastaPart
    {
        private final Path mPath;
        private final BufferedWriter mWriter;
        private int mReadCount;

        static FastaPart create(String outputFastaFile, int index)
        {
            String path = outputFastaFile + ".part" + index;
            try
            {
                return new FastaPart(path, createBufferedWriter(path));
            }
            catch(IOException e)
            {
                throw new RuntimeException("Failed to create candidate FASTA part", e);
            }
        }

        private FastaPart(String path, BufferedWriter writer)
        {
            mPath = new File(path).toPath();
            mWriter = writer;
        }

        private void add(SAMRecord record)
        {
            try
            {
                writeFastaRecord(mWriter, record);
            }
            catch(IOException e)
            {
                throw new RuntimeException("Failed to write candidate reads", e);
            }
            ++mReadCount;
        }

        private Path path() { return mPath; }

        private int readCount() { return mReadCount; }

        private void close() throws IOException
        {
            mWriter.close();
        }
    }

    private static int concatFastaParts(List<FastaPart> parts, String outputFastaFile)
    {
        long startTimeMs = System.currentTimeMillis();
        int candidateCount = 0;
        try(OutputStream out = new BufferedOutputStream(new FileOutputStream(outputFastaFile)))
        {
            for(FastaPart part : parts)
            {
                part.close();
                Files.copy(part.path(), out);
                Files.delete(part.path());
                candidateCount += part.readCount();
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to join candidate FASTA parts", e);
        }

        LOGGER.debug("Joined {} FASTA parts in {}s", parts.size(), secondsSince(startTimeMs));
        return candidateCount;
    }

    private static void closeReaders(List<SamReader> readers)
    {
        for(SamReader reader : readers)
        {
            try
            {
                reader.close();
            }
            catch(IOException e)
            {
                throw new RuntimeException("failed to close tumor BAM", e);
            }
        }
    }
}
