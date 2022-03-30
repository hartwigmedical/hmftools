package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFactory;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountSupplier
{
    private final int mWindowSize;
    private final int mMinMappingQuality;
    private final ExecutorService mExecutorService;
    private final SamReaderFactory mReaderFactory;

    public CountSupplier(
            final int windowSize, final int minMappingQuality,
            final ExecutorService executorService, final SamReaderFactory readerFactory)
    {
        mWindowSize = windowSize;
        mMinMappingQuality = minMappingQuality;
        mExecutorService = executorService;
        mReaderFactory = readerFactory;
    }

    public Multimap<Chromosome, CobaltCount> generateCounts(
            @Nullable final String referenceBam, @Nullable final String tumorBam)
            throws IOException, ExecutionException, InterruptedException
    {
        if (referenceBam == null && tumorBam == null)
        {
            CB_LOGGER.error("No bam file supplied");
            return null;
        }

        final List<ChromosomeLength> lengths;
        try(SamReader reader = mReaderFactory.open(new File(tumorBam != null ? tumorBam : referenceBam)))
        {
            lengths = ChromosomeLengthFactory.create(reader.getFileHeader());
        }

        List<Future<ChromosomeReadCount>> tumorFutures = null;
        List<Future<ChromosomeReadCount>> referenceFutures = null;
        Multimap<Chromosome, ReadCount> tumorCounts = null;
        Multimap<Chromosome, ReadCount> referenceCounts = null;

        if (tumorBam != null)
        {
            CB_LOGGER.info("Calculating Read Count from {}", tumorBam);
            tumorFutures = createFutures(mReaderFactory, new File(tumorBam), lengths);
        }

        if (referenceBam != null)
        {
            CB_LOGGER.info("Calculating Read Count from {}", referenceBam);
            referenceFutures = createFutures(mReaderFactory, new File(referenceBam), lengths);
        }

        if (tumorFutures != null)
        {
            tumorCounts = fromFutures(tumorFutures);
        }
        if (referenceFutures != null)
        {
            referenceCounts = fromFutures(referenceFutures);
        }

        CB_LOGGER.info("Read Count Complete");

        if (tumorBam != null)
        {
            if (referenceBam != null)
            {
                return CobaltCountFactory.pairedTumorNormal(referenceCounts, tumorCounts);
            }
            else
            {
                return CobaltCountFactory.tumorOnly(tumorCounts);
            }
        }
        else
        {
            return CobaltCountFactory.referenceOnly(referenceCounts);
        }
    }

    private List<Future<ChromosomeReadCount>> createFutures(
            final SamReaderFactory readerFactory, final File file, final List<ChromosomeLength> lengths)
    {
        final List<Future<ChromosomeReadCount>> futures = new ArrayList<>();
        for(ChromosomeLength chromosome : lengths)
        {
            final ChromosomeReadCount callable = new ChromosomeReadCount(
                    file, readerFactory, chromosome.chromosome(), mWindowSize, mMinMappingQuality);

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
