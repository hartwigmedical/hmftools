package com.hartwig.hmftools.id;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class HmfIdConfig
{
    public final String Password;
    public final String NewPassword;
    public final String InputHashFile;
    public final String OutputHashFile;
    public final boolean Restore;
    public final int MaxPrecomputeCount;

    public static final Logger ID_LOGGER = LogManager.getLogger(HmfIdConfig.class);

    public static final String PASSWORD = "password";
    public static final String NEW_PASSWORD = "new_password";
    public static final String HASH_FILE_IN = "input_sample_file";
    public static final String HASH_FILE_OUT = "output_sample_file";
    public static final String RESTORE = "restore";
    public static final String MAX_SAMPLE_COUNT = "max_sample_count";

    public static final String DATA_DELIM = ",";
    private static final int DEFAULT_PRECOMPUTE_COUNT = 50000;

    public HmfIdConfig(final ConfigBuilder configBuilder)
    {
        Password = configBuilder.getValue(PASSWORD);
        NewPassword = configBuilder.getValue(NEW_PASSWORD, Password);
        InputHashFile = configBuilder.getValue(HASH_FILE_IN);
        OutputHashFile = configBuilder.getValue(HASH_FILE_OUT);
        Restore = configBuilder.hasFlag(RESTORE);
        MaxPrecomputeCount = configBuilder.getInteger(MAX_SAMPLE_COUNT);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(PASSWORD, true, "Password for ?");
        configBuilder.addConfigItem(NEW_PASSWORD, false, "New password");
        configBuilder.addPath(HASH_FILE_IN, true, "Input sample hash file");
        configBuilder.addConfigItem(HASH_FILE_OUT, true, "Output sample hash file");
        configBuilder.addFlag(RESTORE, "Restore routine");

        configBuilder.addInteger(MAX_SAMPLE_COUNT, "Maximum supported precomputed sample hash count", DEFAULT_PRECOMPUTE_COUNT);

        addDatabaseCmdLineArgs(configBuilder, true);
        ConfigUtils.addLoggingOptions(configBuilder);
    }

}
