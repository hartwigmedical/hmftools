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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Scans the tumor BAM/CRAM once and writes candidate viral reads single-end to a gzipped FASTA, each read once.
// FASTA ids carry the mate number so the two reads of a pair stay distinct. With multiple threads the scan is
// sharded by genome partition (plus the unmapped tail), which requires an indexed input; each worker writes its
// own gzipped part lock-free and the parts are concatenated into the output once every task has finished.
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

        LOGGER.info("extracted {} candidate reads to {}", candidateCount, outputFastaFile);
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
                if(isExcluded(record))
                {
                    continue;
                }

                if(mFilter.isCandidate(record))
                {
                    writeFasta(writer, record);
                    ++candidateCount;
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to extract candidate reads", e);
        }

        return candidateCount;
    }

    private int extractParallel(SamReaderFactory factory, String tumorBamFile, String outputFastaFile)
    {
        SAMSequenceDictionary dictionary;
        try(SamReader reader = factory.open(new File(tumorBamFile)))
        {
            if(!reader.hasIndex())
            {
                throw new UserInputError("multi-threaded extraction requires an indexed BAM/CRAM: " + tumorBamFile);
            }
            dictionary = reader.getFileHeader().getSequenceDictionary();
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to open tumor BAM", e);
        }

        List<ChrBaseRegion> partitions = new ArrayList<>();
        for(SAMSequenceRecord sequence : dictionary.getSequences())
        {
            partitions.addAll(partitionChromosome(sequence, EXTRACTION_PARTITION_SIZE));
        }

        BamSlicer slicer = new BamSlicer(0, true, true, true);
        slicer.setKeepUnmapped();

        // One reader and one gzipped FASTA part per worker thread, both created lazily on first use, so scanning
        // and compression run lock-free. Readers are closed and parts concatenated once every task has finished.
        List<SamReader> readers = Collections.synchronizedList(new ArrayList<>());
        List<FastaPart> parts = Collections.synchronizedList(new ArrayList<>());
        AtomicInteger nextPartIndex = new AtomicInteger();

        ThreadLocal<SamReader> threadReader = ThreadLocal.withInitial(() -> {
            SamReader reader = factory.open(new File(tumorBamFile));
            readers.add(reader);
            return reader;
        });
        ThreadLocal<FastaPart> threadPart = ThreadLocal.withInitial(() -> {
            FastaPart part = FastaPart.create(outputFastaFile, nextPartIndex.getAndIncrement());
            parts.add(part);
            return part;
        });

        ExecutorService executor = Executors.newFixedThreadPool(mThreads);
        List<CompletableFuture<Void>> futures = new ArrayList<>();

        for(ChrBaseRegion region : partitions)
        {
            futures.add(CompletableFuture.runAsync(() -> sliceRegion(slicer, threadReader, threadPart, region), executor));
        }

        // The unmapped reads sit in a single block that cannot be sharded, so they are scanned by one task; its
        // duration relative to the region tasks shows whether the unmapped tail dominates the run.
        futures.add(CompletableFuture.runAsync(() -> sliceUnmapped(slicer, threadReader, threadPart), executor));

        try
        {
            CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
        }
        catch(InterruptedException | ExecutionException e)
        {
            throw new RuntimeException("candidate extraction failed", e);
        }
        finally
        {
            executor.shutdown();
        }

        closeReaders(readers);
        return joinParts(parts, outputFastaFile);
    }

    private void sliceRegion(BamSlicer slicer, ThreadLocal<SamReader> threadReader, ThreadLocal<FastaPart> threadPart, ChrBaseRegion region)
    {
        long startTimeMs = System.currentTimeMillis();
        FastaPart part = threadPart.get();
        int startCount = part.count();

        slicer.slice(threadReader.get(), region, record -> {
            if(isExcluded(record))
            {
                return;
            }

            // A mapped read is owned by the partition containing its start, so copies returned by an overlapping
            // neighbour partition are ignored.
            if(record.getAlignmentStart() < region.start())
            {
                return;
            }

            if(mFilter.isCandidate(record))
            {
                part.add(record);
            }
        });

        LOGGER.debug("region({}) {} candidates in {}s", region, part.count() - startCount, secondsSince(startTimeMs));
    }

    private void sliceUnmapped(BamSlicer slicer, ThreadLocal<SamReader> threadReader, ThreadLocal<FastaPart> threadPart)
    {
        long startTimeMs = System.currentTimeMillis();
        FastaPart part = threadPart.get();
        int startCount = part.count();

        slicer.queryUnmapped(threadReader.get(), record -> {
            if(isExcluded(record))
            {
                return;
            }
            if(mFilter.isCandidate(record))
            {
                part.add(record);
            }
        });

        LOGGER.debug("unmapped reads {} candidates in {}s", part.count() - startCount, secondsSince(startTimeMs));
    }

    private static int joinParts(List<FastaPart> parts, String outputFastaFile)
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
                candidateCount += part.count();
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to join candidate FASTA parts", e);
        }

        LOGGER.debug("joined {} FASTA parts in {}s", parts.size(), secondsSince(startTimeMs));
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

    private static boolean isExcluded(SAMRecord record)
    {
        return record.getDuplicateReadFlag() || record.isSecondaryOrSupplementary();
    }

    private static void writeFasta(BufferedWriter writer, SAMRecord record) throws IOException
    {
        writer.write(">" + fastaId(record));
        writer.newLine();
        writer.write(record.getReadString());
        writer.newLine();
    }

    private static String fastaId(SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            return record.getReadName() + (record.getFirstOfPairFlag() ? "/1" : "/2");
        }
        return record.getReadName();
    }

    private static String secondsSince(long startTimeMs)
    {
        return String.format("%.1f", (System.currentTimeMillis() - startTimeMs) / 1000.0);
    }

    // A worker's private gzipped FASTA shard, written lock-free by the single thread that owns it. Independently
    // valid gzip, so the shards concatenate byte-wise into one gzip stream. Not thread-safe: one instance per thread.
    private static class FastaPart
    {
        private final Path mPath;
        private final BufferedWriter mWriter;
        private int mCount;

        static FastaPart create(String outputFastaFile, int index)
        {
            String base = outputFastaFile.endsWith(".gz") ? outputFastaFile.substring(0, outputFastaFile.length() - 3) : outputFastaFile;
            String path = base + ".part" + index + ".gz";
            try
            {
                return new FastaPart(path, createBufferedWriter(path));
            }
            catch(IOException e)
            {
                throw new RuntimeException("failed to create candidate FASTA part", e);
            }
        }

        private FastaPart(String path, BufferedWriter writer)
        {
            mPath = new File(path).toPath();
            mWriter = writer;
        }

        void add(SAMRecord record)
        {
            try
            {
                writeFasta(mWriter, record);
            }
            catch(IOException e)
            {
                throw new RuntimeException("failed to write candidate reads", e);
            }
            ++mCount;
        }

        void close() throws IOException
        {
            mWriter.close();
        }

        Path path() { return mPath; }

        int count() { return mCount; }
    }
}
