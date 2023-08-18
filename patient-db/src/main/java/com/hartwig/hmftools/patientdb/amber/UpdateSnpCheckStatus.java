package com.hartwig.hmftools.patientdb.amber;

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

public class UpdateSnpCheckStatus
{
    private static final String SAMPLE = "sample";
    private static final String IS_PASS = "is_pass";

    private static final String HAS_PASSED = "true";
    private static final String HAS_FAILED = "false";

    public static void main(@NotNull String[] args) throws ParseException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(IS_PASS, true, "Either pass 'true' or 'false' to set snpcheck pass status");
        addDatabaseCmdLineArgs(configBuilder, true);

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);
        String isPassString = configBuilder.getValue(IS_PASS);

        if(!isPassString.equals(HAS_PASSED) && !isPassString.equals(HAS_FAILED))
        {
            LOGGER.error("pass-string({}) must have value of true or false", isPassString);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(configBuilder);

        boolean isPass = isPassString.equals(HAS_PASSED);
        LOGGER.info("SnpCheck for sample({}) passStatus({})", sample, isPass);
        dbWriter.writeSnpCheck(sample, isPass);

        LOGGER.info("SnpCheck complete");
    }
}
