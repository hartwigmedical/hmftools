package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Multimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.cobalt.count.CountSupplier;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.ReadCountFile;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.gc.GCProfileFactory;
import com.hartwig.hmftools.common.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CountBamLinesApplication {
    private static final Logger LOGGER = LogManager.getLogger(CountBamLinesApplication.class);

    private static final String THREADS = "threads";
    private static final String REFERENCE = "reference";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR = "tumor";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String GC_PROFILE = "gc_profile";
    private static final String WINDOW_SIZE = "window_size";
    private static final String MIN_MAPPING_QUALITY = "min_quality";

    private static final int WINDOW_SIZE_DEFAULT = 1000;
    private static final int MIN_MAPPING_QUALITY_DEFAULT = 10;

    public static void main(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        new CountBamLinesApplication(args);
    }

    private CountBamLinesApplication(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        final VersionInfo versionInfo = new VersionInfo("cobalt.version");
        LOGGER.info("COBALT version: {}", versionInfo.version());

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(GC_PROFILE) || !cmd.hasOption(OUTPUT_DIR) || !cmd.hasOption(TUMOR) || !cmd.hasOption(REFERENCE)) {
            printUsageAndExit(options);
        }

        final String tumor = cmd.getOptionValue(TUMOR);
        final String reference = cmd.getOptionValue(REFERENCE);
        final String outputDirectory = cmd.getOptionValue(OUTPUT_DIR);
        final String outputFilename = CobaltRatioFile.generateFilename(outputDirectory, tumor);
        final File outputFile = new File(outputFilename);

        final String referenceReadCountFilename = ReadCountFile.generateFilename(outputDirectory, reference);
        final File referenceReadCountFile = new File(referenceReadCountFilename);
        final String tumorReadCountFilename = ReadCountFile.generateFilename(outputDirectory, tumor);
        final File tumorReadCountFile = new File(tumorReadCountFilename);

        final boolean existingCobaltOutput = outputFile.exists();
        final boolean existingCobaltReadCounts = referenceReadCountFile.exists() && tumorReadCountFile.exists();
        final boolean useBams = cmd.hasOption(REFERENCE_BAM) || cmd.hasOption(TUMOR_BAM);

        if (useBams) {
            if (!cmd.hasOption(TUMOR_BAM)) {
                LOGGER.warn("Argument -tumor_bam must be set when using -reference_bam");
                printUsageAndExit(options);
            }
            if (!cmd.hasOption(REFERENCE_BAM)) {
                LOGGER.warn("Argument -reference_bam must be set when using -tumor_bam");
                printUsageAndExit(options);
            }
        } else if (!existingCobaltOutput && !existingCobaltReadCounts) {
            LOGGER.warn("BAM arguments may only be omitted when previous cobalt output available.");
            printUsageAndExit(options);
        }

        final Path outputPath = selectOrCreateOutputPath(outputDirectory);
        if (outputPath == null) {
            System.exit(1);
        }

        final int threadCount = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : 4;
        final int windowSize = cmd.hasOption(WINDOW_SIZE) ? Integer.valueOf(cmd.getOptionValue(WINDOW_SIZE)) : WINDOW_SIZE_DEFAULT;
        final int minMappingQuality =
                cmd.hasOption(MIN_MAPPING_QUALITY) ? Integer.valueOf(cmd.getOptionValue(MIN_MAPPING_QUALITY)) : MIN_MAPPING_QUALITY_DEFAULT;
        LOGGER.info("Thread Count: {}, Window Size: {}, Min Quality {}", threadCount, windowSize, minMappingQuality);

        LOGGER.info("Reading GC Profile");
        final Multimap<String, GCProfile> gcProfiles = GCProfileFactory.loadGCContent(windowSize, cmd.getOptionValue(GC_PROFILE));

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("bam-%d").build();
        final ExecutorService executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory);

        final CountSupplier countSupplier =
                new CountSupplier(reference, tumor, outputDirectory, windowSize, minMappingQuality, executorService);
        final Multimap<Chromosome, CobaltCount> readCounts;
        if (useBams) {
            readCounts = countSupplier.fromBam(cmd.getOptionValue(REFERENCE_BAM), cmd.getOptionValue(TUMOR_BAM));
        } else if (existingCobaltOutput) {
            readCounts = countSupplier.fromExistingCobaltFile();
        } else {
            readCounts = countSupplier.fromReadCountFiles();
        }

        final RatioSupplier ratioSupplier = new RatioSupplier(reference, tumor, outputDirectory);
        final Multimap<Chromosome, CobaltRatio> ratios = ratioSupplier.generateRatios(gcProfiles, readCounts);
        LOGGER.info("Persisting cobalt ratios to {}", outputFilename);
        versionInfo.write(outputDirectory);
        CobaltRatioFile.write(outputFilename, ratios);

        new RatioSegmentation(executorService, outputDirectory).applySegmentation(reference, tumor);

        executorService.shutdown();
    }

    @Nullable
    private static Path selectOrCreateOutputPath(@NotNull final String outputString) {
        final File outputFile = new File(outputString);
        if (!outputFile.exists() && !outputFile.mkdirs()) {
            LOGGER.error("Unable to create output directory {} ", outputString);
            return null;
        }

        if (!outputFile.isDirectory()) {
            LOGGER.error("Output dir {} is not a valid directory", outputString);
            return null;
        }

        return outputFile.toPath();
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("COBALT", options);
        System.exit(1);
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        try {
            return parser.parse(options, args);
        } catch (ParseException e) {
            printUsageAndExit(options);
            throw e;
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(WINDOW_SIZE, true, "Window size. Default 1000.");
        options.addOption(THREADS, true, "Number of threads. Default 4.");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample.");
        options.addOption(TUMOR_BAM, true, "Tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(MIN_MAPPING_QUALITY, true, "Min quality. Default 10.");
        options.addOption(GC_PROFILE, true, "Location of GC Profile.");

        return options;
    }
}
