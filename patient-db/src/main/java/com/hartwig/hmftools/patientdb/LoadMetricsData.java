package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.metrics.WGSMetricQC;
import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
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

    private static final String SAMPLE = "sample";

    private static final String REF_METRICS_FILE = "ref_metrics_file";
    private static final String TUMOR_METRICS_FILE = "tumor_metrics_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);

        String sample = cmd.getOptionValue(SAMPLE);
        String refMetricsFile = cmd.getOptionValue(REF_METRICS_FILE);
        String tumorMetricsFile = cmd.getOptionValue(TUMOR_METRICS_FILE);

        if (Utils.anyNull(userName, password, databaseUrl, sample, refMetricsFile, tumorMetricsFile)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load metrics data", options);
        } else {
            String jdbcUrl = "jdbc:" + databaseUrl;
            DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

            LOGGER.info("Extracting and writing metrics for {}", sample);

            WGSMetrics metrics = WGSMetricsFile.read(refMetricsFile, tumorMetricsFile);
            WGSMetricWithQC wgsMetricWithQC = WGSMetricQC.buildWithQCMetric(metrics);
            dbWriter.writeMetrics(sample, wgsMetricWithQC);
        }
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Sample for which we are going to load the metrics");
        options.addOption(REF_METRICS_FILE, true, "Path towards the metrics file holding the ref sample metrics");
        options.addOption(TUMOR_METRICS_FILE, true, "Path towards the metrics file holding the tumor sample metrics");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }
}
