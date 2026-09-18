package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
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
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Extracts potentially viral reads from the tumor BAM and writes them to single-ended FASTA ready for alignment.
// Workers claim regions from a shared queue, each writing its own FASTA shard which are joined at the end.
public class CandidateReadExtractor
{
    @Nullable
    private final String mRefGenomeFile; // required only for CRAM decode
    private final CandidateReadFilter mFilter;
    private final int mThreads;

    private static final Logger LOGGER = LogManager.getLogger(CandidateReadExtractor.class);

    // Queued in place of a mapped region, standing for the BAM's block of unmapped reads.
    private static final ChrBaseRegion UNMAPPED_READS = new ChrBaseRegion("unmapped", NO_POSITION, NO_POSITION);

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

        Queue<ChrBaseRegion> regions = scanRegions(readerFactory, tumorBamFile);
        TaskQueue<ChrBaseRegion> queue = new TaskQueue<>(regions, "regions", 0);
        List<Worker> workers = new ArrayList<>();

        try
        {
            for(int i = 0; i < min(mThreads, regions.size()); ++i)
            {
                workers.add(new Worker(queue, readerFactory.open(new File(tumorBamFile)), FastaPart.create(outputFastaFile, i)));
            }

            if(!executeRunnables(workers, workers.size()))
            {
                throw new RuntimeException("Candidate read extraction failed");
            }

            int candidateCount = joinFastaParts(workers, outputFastaFile);
            LOGGER.info("Extracted {} candidate reads to {}", candidateCount, outputFastaFile);
            return candidateCount;
        }
        finally
        {
            workers.forEach(Worker::close);
        }
    }

    // The unmapped block is one long scan, so it leads the queue and runs alongside the region scans rather than
    // tacking its full duration onto the end.
    private static Queue<ChrBaseRegion> scanRegions(SamReaderFactory readerFactory, String tumorBamFile)
    {
        try(SamReader reader = readerFactory.open(new File(tumorBamFile)))
        {
            if(!reader.hasIndex())
            {
                throw new UserInputError("Tumor BAM/CRAM is not indexed: " + tumorBamFile);
            }

            Queue<ChrBaseRegion> regions = new ConcurrentLinkedQueue<>();
            regions.add(UNMAPPED_READS);
            for(SAMSequenceRecord sequence : reader.getFileHeader().getSequenceDictionary().getSequences())
            {
                regions.addAll(partitionChromosome(sequence, EXTRACTION_PARTITION_SIZE));
            }
            return regions;
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to open tumor BAM", e);
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
                worker.mPart.close();
                Files.copy(worker.mPart.path(), out);
                candidateCount += worker.mPart.readCount();
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
    private class Worker implements Runnable
    {
        private final TaskQueue<ChrBaseRegion> mRegions;
        private final SamReader mReader;
        private final FastaPart mPart;
        private final BamSlicer mSlicer;

        private Worker(TaskQueue<ChrBaseRegion> regions, SamReader reader, FastaPart part)
        {
            mRegions = regions;
            mReader = reader;
            mPart = part;
            mSlicer = new BamSlicer(0, true, true, true);
            mSlicer.setKeepUnmapped();
        }

        @Override
        public void run()
        {
            while(true)
            {
                try
                {
                    scan(mRegions.removeItem());
                }
                catch(NoSuchElementException e)
                {
                    return;
                }
            }
        }

        private void scan(ChrBaseRegion region)
        {
            long startTimeMs = System.currentTimeMillis();
            int startCount = mPart.readCount();

            if(region == UNMAPPED_READS)
            {
                mSlicer.queryUnmapped(mReader, this::addIfCandidate);
            }
            else
            {
                // A mapped read is owned by the partition containing its start, so copies returned by an overlapping
                // neighbour partition are ignored.
                mSlicer.slice(mReader, region, record ->
                {
                    if(record.getAlignmentStart() >= region.start())
                    {
                        addIfCandidate(record);
                    }
                });
            }

            LOGGER.debug("region({}) {} candidates in {}s",
                    region, mPart.readCount() - startCount, format("%.1f", secondsSinceNow(startTimeMs)));
        }

        private void addIfCandidate(SAMRecord record)
        {
            if(!record.getDuplicateReadFlag() && !record.isSecondaryOrSupplementary() && mFilter.isCandidate(record))
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
