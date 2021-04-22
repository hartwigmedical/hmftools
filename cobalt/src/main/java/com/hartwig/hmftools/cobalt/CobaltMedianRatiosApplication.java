package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CobaltMedianRatiosApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(CobaltMedianRatiosApplication.class);

    private static final String TUMOR = "tumor";
    private static final String REFERENCE = "reference";
    private static final String OUTPUT_DIR = "output_dir";

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(OUTPUT_DIR, true, "Output directory");

        return options;
    }

    public static void main(final String... args) throws IOException {
        final Options options = createOptions();
        try (final CobaltMedianRatiosApplication application = new CobaltMedianRatiosApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CountBamLinesApplication", options);
            System.exit(1);
        }
    }

    private final String inputFile;
    private final String outputFile;

    private CobaltMedianRatiosApplication(final Options options, final String... args) throws ParseException {
        VersionInfo versionInfo = new VersionInfo("cobalt.version");
        LOGGER.info("COBALT version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);
        if (!cmd.hasOption(TUMOR)) {
            throw new ParseException("Missing " + TUMOR + " argument");
        }

        if (!cmd.hasOption(REFERENCE)) {
            throw new ParseException("Missing " + REFERENCE + " argument");
        }

        if (!cmd.hasOption(OUTPUT_DIR)) {
            throw new ParseException("Missing " + OUTPUT_DIR + " argument");
        }

        final String tumor = cmd.getOptionValue(TUMOR);
        final String reference = cmd.getOptionValue(REFERENCE);
        final String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        inputFile = CobaltRatioFile.generateFilenameForReading(outputDir, tumor);
        outputFile = MedianRatioFile.generateFilename(outputDir, reference);

        if (!new File(inputFile).exists()) {
            throw new ParseException("Unable to locate file: " + inputFile);
        }

    }

    public void run() throws IOException {
        LOGGER.info("Reading ratio file: {}", inputFile);
        final ListMultimap<Chromosome, CobaltRatio> ratios = CobaltRatioFile.read(inputFile);

        final List<MedianRatio> medianRatios = MedianRatioFactory.create(ratios);
        LOGGER.info("Writing ratio median file: {}", outputFile);
        MedianRatioFile.write(outputFile, medianRatios);
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close() {
        LOGGER.info("Complete");
    }
}
