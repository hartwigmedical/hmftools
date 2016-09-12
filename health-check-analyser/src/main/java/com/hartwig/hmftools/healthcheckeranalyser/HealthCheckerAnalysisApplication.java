package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.BufferedWriter;
import java.io.FileWriter;
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

    private static final String CSV_OUT_ARGS_DESC = "The file to write csv results to.";
    private static final String CSV_OUT = "csvout";

    @NotNull
    private final String reportsPath;
    @NotNull
    private final String csvOut;

    public static void main(final String... args) throws ParseException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String reportsPath = cmd.getOptionValue(REPORTS_PATH);
        String csvOut = cmd.getOptionValue(CSV_OUT);

        if (reportsPath == null || csvOut == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks-Analyser", options);
            System.exit(1);
        }

        new HealthCheckerAnalysisApplication(reportsPath, csvOut);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(REPORTS_PATH, true, REPORTS_PATH_ARGS_DESC);
        options.addOption(CSV_OUT, true, CSV_OUT_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    HealthCheckerAnalysisApplication(@NotNull final String reportsPath, @NotNull final String csvOut) {
        this.reportsPath = reportsPath;
        this.csvOut = csvOut;
    }

    void runAnalysis() throws IOException {
        HealthCheckReport report = HealthCheckReader.readHealthCheckOutput(reportsPath);
        LOGGER.info(HealthCheckDataToCSV.header(report));
        LOGGER.info(HealthCheckDataToCSV.refSample(report));
        LOGGER.info(HealthCheckDataToCSV.tumorSample(report));

        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOut, false));
        writer.write(HealthCheckDataToCSV.header(report));
        writer.newLine();
        writer.write(HealthCheckDataToCSV.refSample(report));
        writer.newLine();
        writer.write(HealthCheckDataToCSV.tumorSample(report));
        writer.close();
    }
}
