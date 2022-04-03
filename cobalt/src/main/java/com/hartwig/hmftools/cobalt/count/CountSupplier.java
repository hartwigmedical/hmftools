package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.Chromosome;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SamReaderFactory;

public class CountSupplier
{
    private final int mWindowSize;
    private final int mMinMappingQuality;
    Multimap<Chromosome, ReadCount> mReferenceCounts = null;
    Multimap<Chromosome, ReadCount> mTumorCounts = null;
    private final ExecutorService mExecutorService;
    private final SamReaderFactory mReaderFactory;
    private final Collection<Chromosome> mChromosomes;

    public Multimap<Chromosome, ReadCount> getReferenceCounts() { return mReferenceCounts; }
    public Multimap<Chromosome, ReadCount> getTumorCounts() { return mTumorCounts; }

    public CountSupplier(
            final int windowSize, final int minMappingQuality,
            final ExecutorService executorService, final SamReaderFactory readerFactory,
            final Collection<Chromosome> chromosomes)
    {
        mWindowSize = windowSize;
        mMinMappingQuality = minMappingQuality;
        mExecutorService = executorService;
        mReaderFactory = readerFactory;
        mChromosomes = chromosomes;
    }

    public void generateCounts(
            @Nullable final String referenceBam, @Nullable final String tumorBam)
            throws ExecutionException, InterruptedException
    {
        if (referenceBam == null && tumorBam == null)
        {
            CB_LOGGER.error("No bam file supplied");
            return;
        }

        List<Future<ChromosomeReadCount>> tumorFutures = null;
        List<Future<ChromosomeReadCount>> referenceFutures = null;

        if (tumorBam != null)
        {
            CB_LOGGER.info("Calculating Read Count from {}", tumorBam);
            tumorFutures = createFutures(mReaderFactory, new File(tumorBam));
        }

        if (referenceBam != null)
        {
            CB_LOGGER.info("Calculating Read Count from {}", referenceBam);
            referenceFutures = createFutures(mReaderFactory, new File(referenceBam));
        }

        if (tumorFutures != null)
        {
            mTumorCounts = fromFutures(tumorFutures);
        }
        if (referenceFutures != null)
        {
            mReferenceCounts = fromFutures(referenceFutures);
        }

        CB_LOGGER.info("Read Count Complete");
    }

    private List<Future<ChromosomeReadCount>> createFutures(
            final SamReaderFactory readerFactory, final File file)
    {
        final List<Future<ChromosomeReadCount>> futures = new ArrayList<>();
        for(Chromosome chromosome : mChromosomes)
        {
            final ChromosomeReadCount callable = new ChromosomeReadCount(
                    file, readerFactory, chromosome, mWindowSize, mMinMappingQuality);

            futures.add(mExecutorService.submit(callable));
        }

        return futures;
    }

    private static Multimap<Chromosome, ReadCount> fromFutures(List<Future<ChromosomeReadCount>> futures)
            throws ExecutionException, InterruptedException
    {
        final ListMultimap<Chromosome, ReadCount> readCounts = ArrayListMultimap.create();
        for(Future<ChromosomeReadCount> future : futures)
        {
            final ChromosomeReadCount readCount = future.get();
            final Chromosome chromosome = readCount.chromosome();
            final List<ReadCount> result = readCount.readCount();

            readCounts.putAll(chromosome, result);
        }

        return readCounts;
    }
}
