package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReportFactory;

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

    public static void main(final String... args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String reportsPath = cmd.getOptionValue(REPORTS_PATH);
        String csvOut = cmd.getOptionValue(CSV_OUT);

        if (reportsPath == null || csvOut == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks-Analyser", options);
            System.exit(1);
        }

        new HealthCheckerAnalysisApplication(reportsPath, csvOut).runAnalysis();
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
        List<HealthCheckReport> reports = Lists.newArrayList();
        File[] runs = new File(reportsPath).listFiles();
        if (runs != null) {
            for (File runDirectory : runs) {
                LOGGER.info("Processing " + runDirectory.getPath());
                HealthCheckReport report = HealthCheckReportFactory.fromHealthCheckReport(runDirectory.getPath());
                reports.add(report);
            }
        } else {
            LOGGER.warn(String.format("%s contains no runs!", reportsPath));
        }

        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOut, false));
        writer.write(HealthCheckDataToCSV.header(reports.get(0)));

        for (HealthCheckReport report : reports) {
            writer.newLine();
            writer.write(HealthCheckDataToCSV.refSample(report));
            writer.newLine();
            writer.write(HealthCheckDataToCSV.tumorSample(report));
        }
        writer.close();
    }
}
