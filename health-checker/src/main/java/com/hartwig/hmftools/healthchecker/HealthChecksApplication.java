package com.hartwig.hmftools.healthchecker;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HealthChecksApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String METRICS_DIR = "metrics_dir";
    private static final String AMBER_DIR = "amber_dir";
    private static final String PURPLE_DIR = "purple_dir";

    private static final String OUTPUT_DIR = "output_dir";

    @NotNull
    private final String refSample;
    @NotNull
    private final String getMetricsDir;

    private HealthChecksApplication(@NotNull final String refSample, @NotNull final String getMetricsDir) {
        this.refSample = refSample;
        this.getMetricsDir = getMetricsDir;
    }

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String refSample = cmd.getOptionValue(REF_SAMPLE);
        final String metricsDir = cmd.getOptionValue(METRICS_DIR);

        if (refSample == null || metricsDir == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks", options);
            System.exit(1);
        }

        HealthChecksApplication app = new HealthChecksApplication(refSample, metricsDir);
        if (cmd.hasOption(TUMOR_SAMPLE)) {
            final String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
            final String amberDir = cmd.getOptionValue(AMBER_DIR);
            final String purpleDir = cmd.getOptionValue(PURPLE_DIR);

            app.runSomatic(tumorSample, amberDir, purpleDir);
        } else {
            app.runSingleSample();
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REF_SAMPLE, true, "The name of the reference sample");
        options.addOption(TUMOR_SAMPLE, true, "The name of the tumor sample");
        options.addOption(PURPLE_DIR, true, "The directory holding the purple output");
        options.addOption(AMBER_DIR, true, "The directory holding the amber output");
        options.addOption(METRICS_DIR, true, "The directory holding the metrics output");

        options.addOption(OUTPUT_DIR, true, "The directory where health checker will write output to");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private void runSingleSample() throws IOException {
        //        final Report report = new JsonReport();
        //        final Collection<HealthChecker> checkers = Lists.newArrayList(new CoverageChecker(), new PurpleChecker(), new AmberChecker());
        //
        //        for (final HealthChecker checker : checkers) {
        //            report.addResult(checker.run(runContext));
        //        }
        //
        //        final Optional<String> reportPath = report.generateReport(reportFilePath);
        //        reportPath.ifPresent(path -> LOGGER.info(String.format("Report generated -> \n%s", path)));
    }

    private void runSomatic(@NotNull String tumorSample, @NotNull String amberDir, @NotNull String purpleDir) throws IOException {
        //        final Report report = new JsonReport();
        //        final Collection<HealthChecker> checkers = Lists.newArrayList(new CoverageChecker(), new PurpleChecker(), new AmberChecker());
        //
        //        for (final HealthChecker checker : checkers) {
        //            report.addResult(checker.run(runContext));
        //        }
        //
        //        final Optional<String> reportPath = report.generateReport(reportFilePath);
        //        reportPath.ifPresent(path -> LOGGER.info(String.format("Report generated -> \n%s", path)));
    }
}
