package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

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
import org.jetbrains.annotations.Nullable;

public class HealthCheckerAnalysisApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckerAnalysisApplication.class);

    private static final String REPORT_PATH_ARGS_DESC = "A path towards a single report to print in STDOUT.";
    private static final String REPORT_PATH = "report";

    private static final String DATA_PATH_ARGS_DESC = "A path which contains multiple health check runs.";
    private static final String DATA_PATH = "data";

    private static final String CSV_OUT_ARGS_DESC = "The file to write csv results to.";
    private static final String CSV_OUT = "csvout";

    @Nullable
    private final String reportPath;
    @Nullable
    private final String dataPath;
    @Nullable
    private final String csvOut;

    public static void main(final String... args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String reportPath = cmd.getOptionValue(REPORT_PATH);
        String dataPath = cmd.getOptionValue(DATA_PATH);
        String csvOut = cmd.getOptionValue(CSV_OUT);

        if (reportPath == null && dataPath == null || (dataPath != null && csvOut == null)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks-Analyser", options);
            System.exit(1);
        }

        new HealthCheckerAnalysisApplication(reportPath, dataPath, csvOut).runAnalysis();
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(REPORT_PATH, true, REPORT_PATH_ARGS_DESC);
        options.addOption(DATA_PATH, true, DATA_PATH_ARGS_DESC);
        options.addOption(CSV_OUT, true, CSV_OUT_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    HealthCheckerAnalysisApplication(@Nullable final String reportPath, @Nullable final String dataPath,
            @Nullable final String csvOut) {
        this.reportPath = reportPath;
        this.dataPath = dataPath;
        this.csvOut = csvOut;
    }

    void runAnalysis() throws IOException {
        if (reportPath == null && dataPath != null && csvOut != null) {
            runInBatchMode(dataPath, csvOut);
        } else if (reportPath != null) {
            runInSingleMode(reportPath);
        } else {
            LOGGER.warn("Could not run, invalid input parameters!");
        }
    }

    private static void runInBatchMode(@NotNull String dataPath, @NotNull String csvOut) throws IOException {
        List<HealthCheckReport> reports = Lists.newArrayList();
        File[] runs = new File(dataPath).listFiles();
        if (runs != null) {
            for (File runDirectory : runs) {
                LOGGER.info("Processing " + runDirectory.getPath());
                HealthCheckReport report = HealthCheckReportFactory.fromHealthCheckReport(runDirectory.getPath());
                reports.add(report);
            }
        } else {
            LOGGER.warn(String.format("%s contains no runs!", dataPath));
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

    private static void runInSingleMode(@NotNull String reportPath) throws FileNotFoundException {
        HealthCheckReport report = HealthCheckReportFactory.fromHealthCheckReport(reportPath);

        LOGGER.info("REF SAMPLE CHECKS - " + report.refSample());
        printChecks(report.refChecks());

        LOGGER.info("");
        LOGGER.info("TUMOR SAMPLE CHECKS - " + report.tumorSample());
        printChecks(report.tumorChecks());

        LOGGER.info("");
        LOGGER.info("PATIENT CHECKS");
        printChecks(report.patientChecks());
    }

    private static void printChecks(@NotNull Map<String, String> checks) {
        for (Map.Entry<String, String> check : checks.entrySet()) {
            LOGGER.info(check.getKey() + "\t" + check.getValue());
        }
    }
}
