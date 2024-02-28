package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.sql.SQLException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CheckSampleAbsent
{
    public static void main(@NotNull String[] args) throws ParseException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        try(DatabaseAccess dbAccess = databaseAccess(configBuilder))
        {
            if(!dbAccess.readPurpleSampleList().contains(sample))
            {
                LOGGER.info("sample({}) is absent from the database", sample);
            }
            else
            {
                LOGGER.warn("sample({}) is not absent from the database", sample);
                System.exit(2);
            }
        }
        catch(Exception e)
        {
            LOGGER.error("failed to check whether sample({}) is absent from the database: {}", sample, e);
            System.exit(1);
        }
    }
}
