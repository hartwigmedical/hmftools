package com.hartwig.hmftools.markdups.umi;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class UmiConfig
{
    public final boolean Enabled;

    private int mUmiLength; // set and accessed in a thread-safe way

    // config options
    private static final String UMI_ENABLED = "umi_enabled";

    public static final char READ_ID_DELIM = ':';
    public static final String READ_ID_DELIM_STR = String.valueOf(READ_ID_DELIM);

    public UmiConfig(boolean enabled)
    {
        Enabled = enabled;
        mUmiLength = 0;
    }

    public static UmiConfig from(final CommandLine cmd)
    {
        return new UmiConfig(
                cmd.hasOption(UMI_ENABLED));
    }

    public String extractUmiId(final String readId)
    {
        int umiLength = checkAndGetUmiLength(readId);
        return extractUmiId(readId, umiLength);
    }

    private synchronized int checkAndGetUmiLength(final String readId)
    {
        if(mUmiLength > 0)
            return mUmiLength;

        String umiId = extractUmiIdFromReadId(readId);
        mUmiLength = umiId.length();
        return mUmiLength;
    }

    public static String extractUmiId(final String readId, final int umiLength)
    {
        return readId.substring(readId.length() - umiLength);
    }

    public static String extractUmiIdFromReadId(final String readId)
    {
        String[] items = readId.split(READ_ID_DELIM_STR, -1);
        return items[items.length - 1];
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(UMI_ENABLED, false, "Use UMIs for duplicates");
    }
}
