package com.hartwig.hmftools.healthchecker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.runners.MetricsChecker;
import com.hartwig.hmftools.healthchecker.runners.PurpleChecker;

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

public final class HealthChecksApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String METRICS_DIR = "metrics_dir";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String OUTPUT_DIR = "output_dir";

    @NotNull
    private final String refSample;
    @Nullable
    private final String tumorSample;
    @NotNull
    private final String metricsDirectory;
    @Nullable
    private final String purpleDirectory;
    @NotNull
    private final String outputDir;

    @VisibleForTesting
    HealthChecksApplication(@NotNull String refSample, @Nullable String tumorSample, @NotNull String metricsDirectory,
            @Nullable String purpleDirectory, @NotNull String outputDir) {
        this.refSample = refSample;
        this.tumorSample = tumorSample;
        this.metricsDirectory = metricsDirectory;
        this.purpleDirectory = purpleDirectory;
        this.outputDir = outputDir;
    }

    public static void main(final String... args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String refSample = cmd.getOptionValue(REF_SAMPLE);
        String metricsDir = cmd.getOptionValue(METRICS_DIR);
        String outputDir = cmd.getOptionValue(OUTPUT_DIR);

        if (refSample == null || metricsDir == null || outputDir == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks", options);
            System.exit(1);
        }

        String tumorSample = cmd.hasOption(TUMOR_SAMPLE) ? cmd.getOptionValue(TUMOR_SAMPLE) : null;
        String purpleDir = cmd.hasOption(PURPLE_DIR) ? cmd.getOptionValue(PURPLE_DIR) : null;

        new HealthChecksApplication(refSample, tumorSample, metricsDir, purpleDir, outputDir).run(true);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(REF_SAMPLE, true, "The name of the reference sample");
        options.addOption(TUMOR_SAMPLE, true, "The name of the tumor sample");
        options.addOption(PURPLE_DIR, true, "The directory holding the purple output");
        options.addOption(METRICS_DIR, true, "The directory holding the metrics output");

        options.addOption(OUTPUT_DIR, true, "The directory where health checker will write output to");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @VisibleForTesting
    void run(boolean writeOutput) throws IOException {
        List<HealthChecker> checkers;
        if (tumorSample == null || purpleDirectory == null) {
            LOGGER.info("Running in SingleSample mode");
            checkers = Lists.newArrayList(new MetricsChecker(refSample, null, metricsDirectory));
        } else {
            LOGGER.info("Running in Somatic mode");
            checkers = Lists.newArrayList(
                    new MetricsChecker(refSample, tumorSample, metricsDirectory),
                    new PurpleChecker(tumorSample, purpleDirectory)
            );
        }

        List<QCValue> qcValues = Lists.newArrayList();
        for (final HealthChecker checker : checkers) {
            qcValues.addAll(checker.run());
        }

        for (QCValue qcValue : qcValues) {
            LOGGER.info("QC '{}' has value '{}'", qcValue.type(), qcValue.value());
        }

        if (HealthCheckEvaluation.isPass(qcValues)) {
            LOGGER.info("Health check evaluation succeeded.");
            if (writeOutput) {
                new FileOutputStream(fileOutputBasePath() + ".HealthCheckSucceeded").close();
            }

        } else {
            LOGGER.info("Health check evaluation failed!");
            if (writeOutput) {
                new FileOutputStream(fileOutputBasePath() + ".HealthCheckFailed").close();
            }
        }
    }

    @NotNull
    private String fileOutputBasePath() {
        String sample = tumorSample != null ? tumorSample : refSample;
        return outputDir + File.separator + sample;
    }
}
