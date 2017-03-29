package com.hartwig.hmftools.healthcheckreader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthcheckreader.model.HealthCheckReport;
import com.hartwig.hmftools.healthcheckreader.model.HealthCheckReportFactory;

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

public class HealthCheckReaderApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckReaderApplication.class);

    private static final String REPORT_PATH = "report";
    private static final String DATA_PATH = "data";
    private static final String CSV_OUT = "csvout";

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String reportPath = cmd.getOptionValue(REPORT_PATH);
        final String dataPath = cmd.getOptionValue(DATA_PATH);
        final String csvOut = cmd.getOptionValue(CSV_OUT);

        final boolean validBatchMode = dataPath != null && csvOut != null;
        final boolean validSingleMode = reportPath != null;

        if (!validBatchMode && !validSingleMode) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Batch-Health-Check-Reader", options);
            System.exit(1);
        }

        new HealthCheckReaderApplication(reportPath, dataPath, csvOut).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(REPORT_PATH, true, "Path towards a single report to print in STDOUT.");
        options.addOption(DATA_PATH, true, "Path which contains multiple health check runs.");
        options.addOption(CSV_OUT, true, "The file to write csv results to in case of batch mode.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Nullable
    private final String reportPath;
    @Nullable
    private final String dataPath;
    @Nullable
    private final String csvOut;

    HealthCheckReaderApplication(@Nullable final String reportPath, @Nullable final String dataPath,
            @Nullable final String csvOut) {
        this.reportPath = reportPath;
        this.dataPath = dataPath;
        this.csvOut = csvOut;
    }

    void run() throws IOException {
        if (reportPath == null && dataPath != null && csvOut != null) {
            runInBatchMode(dataPath, csvOut);
        } else if (reportPath != null) {
            runInSingleMode(reportPath);
        } else {
            LOGGER.warn("Could not run, invalid input parameters!");
        }
    }

    private static void runInBatchMode(@NotNull final String dataPath, @NotNull final String csvOut)
            throws IOException {
        final List<HealthCheckReport> reports = Lists.newArrayList();
        final File[] runs = new File(dataPath).listFiles();
        if (runs != null) {
            for (final File runDirectory : runs) {
                LOGGER.info("Processing " + runDirectory.getPath());
                reports.add(HealthCheckReportFactory.fromHealthCheckReport(runDirectory.getPath()));
            }
        } else {
            LOGGER.warn(String.format("%s contains no runs!", dataPath));
        }

        final BufferedWriter writer = new BufferedWriter(new FileWriter(csvOut, false));
        writer.write(HealthCheckDataToCSV.header(reports.get(0)));

        for (final HealthCheckReport report : reports) {
            writer.newLine();
            writer.write(HealthCheckDataToCSV.refSample(report));
            writer.newLine();
            writer.write(HealthCheckDataToCSV.tumorSample(report));
        }
        writer.close();
        LOGGER.info("Written results to " + csvOut);
    }

    private static void runInSingleMode(@NotNull final String reportPath) throws FileNotFoundException {
        final HealthCheckReport report = HealthCheckReportFactory.fromHealthCheckReport(reportPath);

        LOGGER.info("REF SAMPLE CHECKS - " + report.refSample());
        printChecks(report.refChecks());

        LOGGER.info("");
        LOGGER.info("TUMOR SAMPLE CHECKS - " + report.tumorSample());
        printChecks(report.tumorChecks());

        LOGGER.info("");
        LOGGER.info("PATIENT CHECKS");
        printChecks(report.patientChecks());
    }

    private static void printChecks(@NotNull final Map<String, String> checks) {
        final Set<String> keys = new TreeSet<>(Comparator.<String>naturalOrder());
        keys.addAll(checks.keySet());
        for (final String key : keys) {
            LOGGER.info(key + "\t" + checks.get(key));
        }
    }
}
