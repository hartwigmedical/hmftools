package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadFlagstatData
{
    private static final String REF_FLAGSTAT_FILE = "ref_flagstat_file";
    private static final String TUMOR_FLAGSTAT_FILE = "tumor_flagstat_file";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(REF_FLAGSTAT_FILE, true, "Path towards the flagstat file holding the ref sample flagstats");
        configBuilder.addPath(TUMOR_FLAGSTAT_FILE, true, "Path towards the flagstat file holding the tumor sample flagstats");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        String refFlagstatFile = configBuilder.getValue(REF_FLAGSTAT_FILE);
        String tumorFlagstatFile = configBuilder.getValue(TUMOR_FLAGSTAT_FILE);

        try (DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Extracting and writing flagstats for {}", sample);

            Flagstat refFlagstat = FlagstatFile.read(refFlagstatFile);
            LOGGER.info(" Read reference sample flagstats from {}", refFlagstatFile);
            Flagstat tumorFlagstat = FlagstatFile.read(tumorFlagstatFile);
            LOGGER.info(" Read tumor sample flagstats from {}", tumorFlagstatFile);

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
