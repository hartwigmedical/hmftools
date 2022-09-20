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

    public static final Logger ID_LOGGER = LogManager.getLogger(HmfIdConfig.class);

    public static final String PASSWORD = "password";
    public static final String NEW_PASSWORD = "new_password";
    public static final String HASH_FILE_IN = "in";
    public static final String HASH_FILE_OUT = "out";
    public static final String RESTORE = "restore";

    public static final String DATA_DELIM = ",";

    public HmfIdConfig(final CommandLine cmd)
    {
        Password = cmd.getOptionValue(PASSWORD);
        NewPassword = cmd.getOptionValue(NEW_PASSWORD);
        InputHashFile = cmd.getOptionValue(HASH_FILE_IN);
        OutputHashFile = cmd.getOptionValue(HASH_FILE_OUT);
        Restore = cmd.hasOption(RESTORE);
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(PASSWORD, true, "Password for ?");
        options.addOption(NEW_PASSWORD, true, "");
        options.addOption(HASH_FILE_IN, true, "");
        options.addOption(HASH_FILE_OUT, true, "");
        options.addOption(RESTORE, false, "");

        addDatabaseCmdLineArgs(options);

        addLoggingOptions(options);
    }

}
