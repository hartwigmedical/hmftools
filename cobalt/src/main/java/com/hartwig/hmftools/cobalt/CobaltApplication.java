package com.hartwig.hmftools.cobalt;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Multimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.count.CountSupplier;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public class CobaltApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(CobaltApplication.class);

    private final CobaltConfig config;
    private final VersionInfo versionInfo;
    private final ExecutorService executorService;

    public static void main(final String... args) throws IOException, ExecutionException, InterruptedException {
        final Options options = CobaltConfig.createOptions();
        try (final CobaltApplication application = new CobaltApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CountBamLinesApplication", options);
            System.exit(1);
        }
    }

    private CobaltApplication(final Options options, final String... args) throws ParseException, IOException {
        versionInfo = new VersionInfo("cobalt.version");
        LOGGER.info("COBALT version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);
        config = CobaltConfig.createConfig(cmd);

        if (!config.refGenomePath().isEmpty() && !new File(config.gcProfilePath()).exists()) {
            throw new IOException("Unable to locate ref genome file " + config.refGenomePath());
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        executorService = Executors.newFixedThreadPool(config.threadCount(), namedThreadFactory);

        LOGGER.info("Thread Count: {}, Window Size: {}, Min Quality {}",
                config.threadCount(),
                config.windowSize(),
                config.minMappingQuality());
    }

    private void run() throws IOException, ExecutionException, InterruptedException {
        LOGGER.info("Reading GC Profile from {}", config.gcProfilePath());
        final Multimap<Chromosome, GCProfile> gcProfiles = GCProfileFactory.loadGCContent(config.windowSize(), config.gcProfilePath());

        final List<BEDFeature> diploidBedFile = diploidBedFile(config);

        final SamReaderFactory readerFactory = readerFactory(config);
        final CountSupplier countSupplier = new CountSupplier(config.tumor(),
                config.outputDirectory(),
                config.windowSize(),
                config.minMappingQuality(),
                executorService,
                readerFactory);
        final Multimap<Chromosome, CobaltCount> readCounts = config.tumorOnly()
                ? countSupplier.tumorOnly(config.tumorBamPath())
                : countSupplier.pairedTumorNormal(config.referenceBamPath(), config.tumorBamPath());

        final RatioSupplier ratioSupplier = new RatioSupplier(config.reference(), config.tumor(), config.outputDirectory());
        final Multimap<Chromosome, CobaltRatio> ratios = config.tumorOnly()
                ? ratioSupplier.tumorOnly(diploidBedFile, gcProfiles, readCounts)
                : ratioSupplier.tumorNormalPair(gcProfiles, readCounts);

        final String outputFilename = CobaltRatioFile.generateFilenameForWriting(config.outputDirectory(), config.tumor());
        LOGGER.info("Persisting cobalt ratios to {}", outputFilename);
        versionInfo.write(config.outputDirectory());
        CobaltRatioFile.write(outputFilename, ratios);

        new RatioSegmentation(executorService, config.outputDirectory()).applySegmentation(config.reference(), config.tumor());
    }

    @NotNull
    private static List<BEDFeature> diploidBedFile(CobaltConfig config) throws IOException {
        List<BEDFeature> result = Lists.newArrayList();
        if (!config.tumorOnly()) {
            return result;
        }

        LOGGER.info("Reading diploid regions from {}", config.tumorOnlyDiploidBed());
        try (final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(config.tumorOnlyDiploidBed(),
                new BEDCodec(),
                false)) {
            for (BEDFeature bedFeature : reader.iterator()) {
                result.add(bedFeature);
            }
        }

        return result;
    }

    @NotNull
    private static SamReaderFactory readerFactory(@NotNull final CobaltConfig config) {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.validationStringency());
        if (!config.refGenomePath().isEmpty()) {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.refGenomePath())));
        }
        return readerFactory;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close() {
        executorService.shutdown();
        LOGGER.info("Complete");
    }
}
