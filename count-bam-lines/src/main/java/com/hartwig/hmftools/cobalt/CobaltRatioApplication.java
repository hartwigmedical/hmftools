package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.concurrent.ExecutionException;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.ratio.RatioSupplier;
import com.hartwig.hmftools.cobalt.segment.PCFSegment;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCountFile;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.gc.GCProfileFactory;

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

public class CobaltRatioApplication {
    private static final Logger LOGGER = LogManager.getLogger(CobaltRatioApplication.class);

    private static final String OUTPUT_DIR = "output_dir";
    private static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String GC_PROFILE = "gc_profile";


    public static void main(final String... args)
            throws ParseException, IOException, ExecutionException, InterruptedException, HartwigException {
        new CobaltRatioApplication(args);
    }

    private CobaltRatioApplication(final String... args)
            throws ParseException, IOException, ExecutionException, InterruptedException, HartwigException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(OUTPUT_DIR) || !cmd.hasOption(SAMPLE)) {
            printUsageAndExit(options);
        }

        final String sample = cmd.getOptionValue(SAMPLE);
        final String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        final Path outputPath = selectOrCreateOutputPath(cmd.getOptionValue(OUTPUT_DIR));
        if (outputPath == null) {
            System.exit(1);
        }

        // Load Data
        LOGGER.info("Reading GC Profile");
        final Multimap<String, GCProfile> gcProfiles = GCProfileFactory.loadGCContent(cmd.getOptionValue(GC_PROFILE));
        final String readCountFilename = ReadCountFile.generateFilename(outputDir, sample);
        LOGGER.info("Reading read count from {}", readCountFilename);
        final Multimap<String, ReadCount>  readCounts = ReadCountFile.readFile(readCountFilename);

        // Generate ratios
        final RatioSupplier ratioSupplier = new RatioSupplier(sample, outputDir, cmd.hasOption(REFERENCE));
        ratioSupplier.generateRatios(gcProfiles, readCounts);

        // Segmentation
        new PCFSegment(outputDir).ratioSegmentation(sample);

    }

    @Nullable
    private Path selectOrCreateOutputPath(@NotNull final String outputString) {
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
        options.addOption(REFERENCE, false, "Apply diploid normalization to ratios");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(GC_PROFILE, true, "Location of GC Profile.");


        return options;
    }
}
