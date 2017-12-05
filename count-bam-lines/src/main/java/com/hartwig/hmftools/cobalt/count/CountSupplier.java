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
import com.hartwig.hmftools.cobalt.CountBamLinesApplication;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFile;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltCountFactory;
import com.hartwig.hmftools.common.cobalt.CobaltCountFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCountFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountSupplier {

    private static final Logger LOGGER = LogManager.getLogger(CountBamLinesApplication.class);

    private final String reference;
    private final String tumor;
    private final String outputDirectory;
    private final int windowSize;
    private final int minMappingQuality;
    private final ExecutorService executorService;

    public CountSupplier(final String reference, final String tumor, final String outputDirectory, final int windowSize,
            final int minMappingQuality, final ExecutorService executorService) {
        this.reference = reference;
        this.tumor = tumor;
        this.outputDirectory = outputDirectory;
        this.windowSize = windowSize;
        this.minMappingQuality = minMappingQuality;
        this.executorService = executorService;
    }

    @NotNull
    public Multimap<Chromosome, CobaltCount> fromReadCountFiles() throws IOException {
        final String referenceFilename = ReadCountFile.generateFilename(outputDirectory, reference);
        LOGGER.info("Reading reference count from {}", referenceFilename);
        Multimap<String, ReadCount> referenceCounts = ReadCountFile.readFile(referenceFilename);

        final String tumorFilename = ReadCountFile.generateFilename(outputDirectory, tumor);
        LOGGER.info("Reading tumor count from {}", tumorFilename);
        Multimap<String, ReadCount> tumorCounts = ReadCountFile.readFile(tumorFilename);

        return CobaltCountFactory.merge(referenceCounts, tumorCounts);
    }

    @NotNull
    public Multimap<Chromosome, CobaltCount> fromExistingCobaltFile() throws IOException {
        final String filename = CobaltRatioFile.generateFilename(outputDirectory, tumor);
        LOGGER.info("Reading reference and tumor counts from {}", filename);
        return CobaltCountFile.read(filename);
    }

    @NotNull
    public Multimap<Chromosome, CobaltCount> fromBam(@NotNull final String referenceBam, @NotNull final String tumorBam)
            throws IOException, ExecutionException, InterruptedException {
        final File tumorFile = new File(tumorBam);
        final File referenceFile = new File(referenceBam);

        final SamReaderFactory readerFactory = SamReaderFactory.make();
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

        final Multimap<String, ReadCount> tumorCounts = fromFutures(tumorFutures);
        final Multimap<String, ReadCount> referenceCounts = fromFutures(referenceFutures);

        LOGGER.info("Read Count Complete");
        return CobaltCountFactory.merge(referenceCounts, tumorCounts);
    }

    @NotNull
    private List<Future<ChromosomeReadCount>> createFutures(final SamReaderFactory readerFactory, final File file, final List<ChromosomeLength> lengths) {
        final List<Future<ChromosomeReadCount>> futures = Lists.newArrayList();
        for (ChromosomeLength chromosome : lengths) {
            final ChromosomeReadCount callable =
                    new ChromosomeReadCount(file, readerFactory, chromosome.chromosome(), chromosome.position(), windowSize,
                            minMappingQuality);
            futures.add(executorService.submit(callable));
        }

        return futures;
    }

    @NotNull
    private static Multimap<String, ReadCount> fromFutures(List<Future<ChromosomeReadCount>> futures)
            throws ExecutionException, InterruptedException {
        final ListMultimap<String, ReadCount> readCounts = ArrayListMultimap.create();
        for (Future<ChromosomeReadCount> future : futures) {
            final ChromosomeReadCount readCount = future.get();
            final String chromosome = readCount.chromosome();
            final List<ReadCount> result = readCount.readCount();

            readCounts.putAll(chromosome, result);
        }

        return readCounts;
    }
}
