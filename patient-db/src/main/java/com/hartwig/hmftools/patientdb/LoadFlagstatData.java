package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadFlagstatData
{
    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, true, REF_METRICS_DIR_DESC);
        configBuilder.addPath(TUMOR_METRICS_DIR_CFG, true, TUMOR_METRICS_DIR_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        String tumor = configBuilder.getValue(SAMPLE);
        String tumorMetricsDir = configBuilder.getValue(TUMOR_METRICS_DIR_CFG);
        String tumorFlagstatFile = BamFlagStats.generateFilename(tumorMetricsDir, tumor);

        String reference = configBuilder.getValue(REFERENCE);
        String refMetricsDir = configBuilder.getValue(REF_METRICS_DIR_CFG);
        String refFlagstatFile = BamFlagStats.generateFilename(refMetricsDir, reference);

        try(DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("load reference sample flagstats from {}", refFlagstatFile);
            BamFlagStats refFlagstat = BamFlagStats.read(refFlagstatFile);

            LOGGER.info("load tumor sample flagstats from {}", tumorFlagstatFile);
            BamFlagStats tumorFlagstat = BamFlagStats.read(tumorFlagstatFile);

            dbWriter.writeFlagstats(sample, refFlagstat, tumorFlagstat);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load flagstats", e);
            System.exit(1);
        }
    }
}
