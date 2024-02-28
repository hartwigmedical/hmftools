package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.metrics.WGSMetricQC;
import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadMetricsData
{
    private static final String REF_METRICS_FILE = "ref_metrics_file";
    private static final String TUMOR_METRICS_FILE = "tumor_metrics_file";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(REF_METRICS_FILE, true, "Path towards the metrics file holding the ref sample metrics");
        configBuilder.addPath(TUMOR_METRICS_FILE, true, "Path towards the metrics file holding the tumor sample metrics");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        String refMetricsFile = configBuilder.getValue(REF_METRICS_FILE);
        String tumorMetricsFile = configBuilder.getValue(TUMOR_METRICS_FILE);

        try (DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Extracting and writing metrics for {}", sample);

            WGSMetrics refMetrics = WGSMetricsFile.read(refMetricsFile);
            LOGGER.info(" Read reference sample metrics from {}", refMetricsFile);
            WGSMetrics tumorMetrics = WGSMetricsFile.read(tumorMetricsFile);
            LOGGER.info(" Read tumor sample metrics from {}", tumorMetricsFile);

            WGSMetricWithQC wgsMetricWithQC = WGSMetricQC.buildWithQCMetric(refMetrics, tumorMetrics);
            dbWriter.writeMetrics(sample, wgsMetricWithQC);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load metrics", e);
            System.exit(1);
        }
    }
}
