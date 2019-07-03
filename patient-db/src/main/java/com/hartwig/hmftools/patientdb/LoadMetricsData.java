package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LoadMetricsData {

    private static final Logger LOGGER = LogManager.getLogger(LoadCanonicalTranscripts.class);

    private static final String METRICS_DIR = "metrics_dir";
    private static final String REF_SAMPLE = "ref_sample";
    private static final String TUMOR_SAMPLE = "tumor_sample";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);

        String metricsDirectory = cmd.getOptionValue(METRICS_DIR);
        String refSample = cmd.getOptionValue(REF_SAMPLE);
        String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);

        if (Utils.anyNull(userName, password, databaseUrl, metricsDirectory, refSample, tumorSample)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load metrics data", options);
        } else {
            final File metricsDir = new File(metricsDirectory);
            if (metricsDir.isDirectory()) {
                final String jdbcUrl = "jdbc:" + databaseUrl;
                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

                LOGGER.info(String.format("Extracting and writing metrics for %s", tumorSample));
                String refFile = WGSMetricsFile.generateFilename(metricsDirectory, refSample);
                String tumorFile = WGSMetricsFile.generateFilename(metricsDirectory, tumorSample);

                WGSMetrics metrics = WGSMetricsFile.read(refFile, tumorFile);
                dbWriter.writeMetrics(tumorSample, metrics);
            } else {
                if (!metricsDir.exists()) {
                    LOGGER.warn("Dir " + metricsDir + " does not exist.");
                }
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("patient-db - load metrics data", options);
            }
        }
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(METRICS_DIR, true, "Path towards the folder containing the metrics files.");
        options.addOption(REF_SAMPLE, true, "ID of the reference sample of set.");
        options.addOption(TUMOR_SAMPLE, true, "ID of the tumor sample of set.");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }
}
