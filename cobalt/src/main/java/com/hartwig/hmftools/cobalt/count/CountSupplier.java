package com.hartwig.hmftools.cobalt.count;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.CobaltApplication;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltCountFactory;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountSupplier {

    private static final Logger LOGGER = LogManager.getLogger(CobaltApplication.class);

    private final String tumor;
    private final String outputDirectory;
    private final int windowSize;
    private final int minMappingQuality;
    private final ExecutorService executorService;
    private final SamReaderFactory readerFactory;

    public CountSupplier(final String tumor, final String outputDirectory, final int windowSize, final int minMappingQuality,
            final ExecutorService executorService, final SamReaderFactory readerFactory) {
        this.tumor = tumor;
        this.outputDirectory = outputDirectory;
        this.windowSize = windowSize;
        this.minMappingQuality = minMappingQuality;
        this.executorService = executorService;
        this.readerFactory = readerFactory;
    }

    @NotNull
    public Multimap<Chromosome, CobaltCount> pairedTumorNormal(@NotNull final String referenceBam, @NotNull final String tumorBam)
            throws IOException, ExecutionException, InterruptedException {
        final File tumorFile = new File(tumorBam);
        final File referenceFile = new File(referenceBam);

        final String chromosomeLengthFileName = ChromosomeLengthFile.generateFilename(outputDirectory, tumor);
        final List<ChromosomeLength> lengths;
        try (SamReader reader = readerFactory.open(tumorFile)) {
            lengths = ChromosomeLengthFactory.create(reader.getFileHeader());
        }
        ChromosomeLengthFile.write(chromosomeLengthFileName, lengths);

        LOGGER.info("Calculating Read Count from {}", tumorFile.toString());
        final List<Future<ChromosomeReadCount>> tumorFutures = createFutures(readerFactory, tumorFile, lengths);

        LOGGER.info("Calculating Read Count from {}", referenceFile.toString());
        final List<Future<ChromosomeReadCount>> referenceFutures = createFutures(readerFactory, referenceFile, lengths);

        final Multimap<Chromosome, ReadCount> tumorCounts = fromFutures(tumorFutures);
        final Multimap<Chromosome, ReadCount> referenceCounts = fromFutures(referenceFutures);

        LOGGER.info("Read Count Complete");
        return CobaltCountFactory.pairedTumorNormal(referenceCounts, tumorCounts);
    }

    @NotNull
    public Multimap<Chromosome, CobaltCount> tumorOnly(@NotNull final String tumorBam)
            throws IOException, ExecutionException, InterruptedException {
        final File tumorFile = new File(tumorBam);

        final String chromosomeLengthFileName = ChromosomeLengthFile.generateFilename(outputDirectory, tumor);
        final List<ChromosomeLength> lengths;
        try (SamReader reader = readerFactory.open(tumorFile)) {
            lengths = ChromosomeLengthFactory.create(reader.getFileHeader());
        }
        ChromosomeLengthFile.write(chromosomeLengthFileName, lengths);

        LOGGER.info("Calculating Read Count from {}", tumorFile.toString());
        final List<Future<ChromosomeReadCount>> tumorFutures = createFutures(readerFactory, tumorFile, lengths);

        final Multimap<Chromosome, ReadCount> tumorCounts = fromFutures(tumorFutures);

        LOGGER.info("Read Count Complete");
        return CobaltCountFactory.tumorOnly(tumorCounts);
    }

    @NotNull
    private List<Future<ChromosomeReadCount>> createFutures(final SamReaderFactory readerFactory, final File file,
            final List<ChromosomeLength> lengths) {
        final List<Future<ChromosomeReadCount>> futures = Lists.newArrayList();
        for (ChromosomeLength chromosome : lengths) {
            final ChromosomeReadCount callable = new ChromosomeReadCount(file,
                    readerFactory,
                    chromosome.chromosome(),
                    chromosome.length(),
                    windowSize,
                    minMappingQuality);
            futures.add(executorService.submit(callable));
        }

        return futures;
    }

    @NotNull
    private static Multimap<Chromosome, ReadCount> fromFutures(List<Future<ChromosomeReadCount>> futures)
            throws ExecutionException, InterruptedException {
        final ListMultimap<Chromosome, ReadCount> readCounts = ArrayListMultimap.create();
        for (Future<ChromosomeReadCount> future : futures) {
            final ChromosomeReadCount readCount = future.get();
            final Chromosome chromosome = readCount.chromosome();
            final List<ReadCount> result = readCount.readCount();

            readCounts.putAll(chromosome, result);
        }

        return readCounts;
    }
}
