package com.hartwig.hmftools.id;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
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

    public HmfIdConfig(final CommandLine cmd)
    {
        Password = cmd.getOptionValue(PASSWORD);
        NewPassword = cmd.getOptionValue(NEW_PASSWORD, Password);
        InputHashFile = cmd.getOptionValue(HASH_FILE_IN);
        OutputHashFile = cmd.getOptionValue(HASH_FILE_OUT);
        Restore = cmd.hasOption(RESTORE);
        MaxPrecomputeCount = Integer.parseInt(cmd.getOptionValue(MAX_SAMPLE_COUNT, String.valueOf(DEFAULT_PRECOMPUTE_COUNT)));
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(PASSWORD, true, "Password for ?");
        options.addOption(NEW_PASSWORD, true, "New password");
        options.addOption(HASH_FILE_IN, true, "Input sample hash file");
        options.addOption(HASH_FILE_OUT, true, "Output sample hash file");
        options.addOption(RESTORE, false, "");

        options.addOption(
                MAX_SAMPLE_COUNT, true,
                "Maximum supported precomputed sample hash count, default: " + DEFAULT_PRECOMPUTE_COUNT);

        addDatabaseCmdLineArgs(options);
        addLoggingOptions(options);
    }

}
