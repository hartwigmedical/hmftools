package com.hartwig.hmftools.common.utils;

import org.apache.commons.cli.CommandLine;

public class ConfigUtils
{
    public static double getConfigValue(final CommandLine cmd, final String configName, double defaultValue)
    {
        return cmd.hasOption(configName) ?  Double.parseDouble(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static int getConfigValue(final CommandLine cmd, final String configName, int defaultValue)
    {
        return cmd.hasOption(configName) ?  Integer.parseInt(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static boolean getConfigValue(final CommandLine cmd, final String configName, boolean defaultValue)
    {
        return cmd.hasOption(configName) ?  Boolean.parseBoolean(cmd.getOptionValue(configName)) : defaultValue;
    }
}
