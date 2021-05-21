package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

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

public class LoadMetricsData {

    private static final Logger LOGGER = LogManager.getLogger(LoadMetricsData.class);

    private static final String SAMPLE = "sample";

    private static final String REF_METRICS_FILE = "ref_metrics_file";
    private static final String TUMOR_METRICS_FILE = "tumor_metrics_file";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String refMetricsFile = cmd.getOptionValue(REF_METRICS_FILE);
        String tumorMetricsFile = cmd.getOptionValue(TUMOR_METRICS_FILE);

        if (Utils.anyNull(sample, refMetricsFile, tumorMetricsFile)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load Metrics Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Extracting and writing metrics for {}", sample);

        WGSMetrics refMetrics = WGSMetricsFile.read(refMetricsFile);
        LOGGER.info(" Read reference sample metrics from {}", refMetricsFile);
        WGSMetrics tumorMetrics = WGSMetricsFile.read(tumorMetricsFile);
        LOGGER.info(" Read tumor sample metrics from {}", tumorMetricsFile);

        WGSMetricWithQC wgsMetricWithQC = WGSMetricQC.buildWithQCMetric(refMetrics, tumorMetrics);
        dbWriter.writeMetrics(sample, wgsMetricWithQC);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Sample for which we are going to load the metrics");
        options.addOption(REF_METRICS_FILE, true, "Path towards the metrics file holding the ref sample metrics");
        options.addOption(TUMOR_METRICS_FILE, true, "Path towards the metrics file holding the tumor sample metrics");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
