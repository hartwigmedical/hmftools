package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.virusdetect.VirusConstants.EXTRACTION_PARTITION_SIZE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.BiConsumer;

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

// Scans the tumor BAM/CRAM once and writes candidate viral reads single-end to a FASTA, each read once.
// FASTA ids carry the mate number so the two reads of a pair stay distinct. With multiple threads the scan
// is sharded by genome partition (plus the unmapped tail), which requires an indexed input.
public class CandidateReadExtractor
{
    @Nullable
    private final String mRefGenomeFile; // required only for CRAM decode
    private final CandidateReadFilter mFilter;
    private final int mThreads;

    private static final Logger LOGGER = LogManager.getLogger(CandidateReadExtractor.class);

    // Per-thread FASTA accumulation flushed to the shared writer once it reaches this many characters.
    private static final int FASTA_BLOCK_FLUSH_CHARS = 256 * 1024;

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

        // Each worker thread formats candidates into its own buffer and takes the shared writer lock only to
        // flush a full block, so lock contention is amortised over many reads. Leftover buffers are flushed
        // single-threaded once every scan task has finished.
        List<FastaBuffer> threadBuffers = Collections.synchronizedList(new ArrayList<>());
        ExecutorService executor = Executors.newFixedThreadPool(mThreads);
        try(BufferedWriter writer = createBufferedWriter(outputFastaFile))
        {
            ThreadLocal<FastaBuffer> threadBuffer = ThreadLocal.withInitial(() -> {
                FastaBuffer buffer = new FastaBuffer(writer);
                threadBuffers.add(buffer);
                return buffer;
            });

            BiConsumer<SAMRecord, ChrBaseRegion> consumer = (record, region) -> {
                if(isExcluded(record))
                {
                    return;
                }

                // A mapped read is owned by the partition containing its start, so copies returned by an
                // overlapping neighbour partition are ignored. Unmapped reads (region == null) are always kept.
                if(region != null && record.getAlignmentStart() < region.start())
                {
                    return;
                }

                if(mFilter.isCandidate(record))
                {
                    threadBuffer.get().add(record);
                }
            };

            slicer.queryAsync(new File(tumorBamFile), factory, partitions, true, executor, consumer).get();

            int candidateCount = 0;
            for(FastaBuffer buffer : threadBuffers)
            {
                candidateCount += buffer.flush();
            }
            return candidateCount;
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to extract candidate reads", e);
        }
        catch(InterruptedException | ExecutionException e)
        {
            throw new RuntimeException("candidate extraction failed", e);
        }
        finally
        {
            executor.shutdown();
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

    // A worker's private accumulation of formatted FASTA records; the shared writer lock is taken only to
    // flush a full block, not per read. Not thread-safe: one instance per thread.
    private static class FastaBuffer
    {
        private final BufferedWriter mWriter;
        private final StringBuilder mPending = new StringBuilder(FASTA_BLOCK_FLUSH_CHARS + 1024);
        private int mCount;

        FastaBuffer(BufferedWriter writer)
        {
            mWriter = writer;
        }

        void add(SAMRecord record)
        {
            mPending.append('>').append(fastaId(record)).append('\n').append(record.getReadString()).append('\n');
            ++mCount;
            if(mPending.length() >= FASTA_BLOCK_FLUSH_CHARS)
            {
                writeBlock();
            }
        }

        int flush()
        {
            writeBlock();
            return mCount;
        }

        private void writeBlock()
        {
            if(mPending.length() == 0)
            {
                return;
            }
            try
            {
                synchronized(mWriter)
                {
                    mWriter.write(mPending.toString());
                }
            }
            catch(IOException e)
            {
                throw new RuntimeException("failed to write candidate reads", e);
            }
            mPending.setLength(0);
        }
    }
}
