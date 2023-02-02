package com.hartwig.hmftools.markdups.umi;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.umi.UmiGroup.exceedsUmiIdDiff;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class UmiConfig
{
    public final boolean Enabled;
    public final boolean Debug;

    private int mUmiLength; // set and accessed in a thread-safe way

    private final Set<String> mDefinedUmis;

    // config options
    private static final String UMI_ENABLED = "umi_enabled";
    private static final String UMI_DEFINED_IDS = "umi_defined_ids";
    private static final String UMI_DEBUG = "umi_debug";

    public static final char READ_ID_DELIM = ':';
    public static final String READ_ID_DELIM_STR = String.valueOf(READ_ID_DELIM);

    public UmiConfig(boolean enabled) { this(enabled, false); }

    public UmiConfig(boolean enabled, boolean debug)
    {
        Enabled = enabled;
        Debug = debug;
        mUmiLength = 0;
        mDefinedUmis = Sets.newHashSet();
    }

    public boolean hasDefinedUmis() { return !mDefinedUmis.isEmpty(); }
    public void addDefinedUmis(final Set<String> umis) { mDefinedUmis.addAll(umis); }

    public static UmiConfig from(final CommandLine cmd)
    {
        UmiConfig umiConfig = new UmiConfig(cmd.hasOption(UMI_ENABLED), cmd.hasOption(UMI_DEBUG));

        String definedUmiIdsFilename = cmd.getOptionValue(UMI_DEFINED_IDS);

        if(definedUmiIdsFilename != null)
        {
            try
            {
                List<String> umis = Files.readAllLines(Paths.get(definedUmiIdsFilename));
                MD_LOGGER.info("loaded {} defined UMI IDs from {}", umis.size(), definedUmiIdsFilename);
                umiConfig.addDefinedUmis(umis.stream().collect(Collectors.toSet()));
            }
            catch(IOException e)
            {
                MD_LOGGER.error("failed to load defined UMI IDs: {}", e.toString());
            }
        }

        return umiConfig;
    }

    public String extractUmiId(final String readId)
    {
        int umiLength = checkAndGetUmiLength(readId);
        return extractUmiId(readId, umiLength);
    }

    public String matchDefinedUmiId(final String umiId)
    {
        if(mDefinedUmis.contains(umiId))
            return umiId;

        for(String definedId : mDefinedUmis)
        {
            if(!exceedsUmiIdDiff(umiId, definedId))
                return definedId;
        }

        return null;
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
        options.addOption(UMI_DEBUG, false, "Debug options for UMIs");
        options.addOption(UMI_DEFINED_IDS, true, "Optional set of defined UMI IDs in file");
    }
}
