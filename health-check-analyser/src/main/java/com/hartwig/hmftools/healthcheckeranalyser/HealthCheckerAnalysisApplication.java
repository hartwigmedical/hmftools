package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.IOException;

import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HealthCheckerAnalysisApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckerAnalysisApplication.class);

    private static final String REPORTS_PATH_ARGS_DESC = "The path where health check report can be found.";
    private static final String REPORTS_PATH = "reports";

    @NotNull
    private final String reportsPath;

    public static void main(final String... args) throws ParseException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String reportsPath = cmd.getOptionValue(REPORTS_PATH);

        if (reportsPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks-Analyser", options);
            System.exit(1);
        }

        new HealthCheckerAnalysisApplication(reportsPath);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(REPORTS_PATH, true, REPORTS_PATH_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    HealthCheckerAnalysisApplication(@NotNull final String reportsPath) {
        this.reportsPath = reportsPath;
    }

    void runAnalysis() throws IOException {
        HealthCheckReport report = HealthCheckReader.readHealthCheckOutput(reportsPath);
        LOGGER.info(HealthCheckDataToCSV.header(report));
        LOGGER.info(HealthCheckDataToCSV.refSample(report));
        LOGGER.info(HealthCheckDataToCSV.tumorSample(report));
    }
}
