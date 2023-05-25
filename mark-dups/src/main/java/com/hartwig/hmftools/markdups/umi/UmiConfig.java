package com.hartwig.hmftools.markdups.umi;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_MAX_UMI_BASE_DIFF;
import static com.hartwig.hmftools.markdups.umi.UmiUtils.exceedsUmiIdDiff;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.markdups.common.Constants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class UmiConfig
{
    public final boolean Enabled;
    public final boolean Duplex; // collapse duplex UMI groups
    public final String DuplexDelim;
    public final boolean CollapseReversed;
    public final boolean BaseStats;
    public final boolean HighlightConsensus; // purely for viewing in IGV
    public final int PermittedBaseDiff;

    private int mUmiLength; // set and accessed in a thread-safe way

    private final Set<String> mDefinedUmis;

    // config options
    private static final String UMI_ENABLED = "umi_enabled";
    private static final String UMI_DUPLEX = "umi_duplex";
    private static final String UMI_DUPLEX_DELIM = "umi_duplex_delim";
    private static final String UMI_DEFINED_IDS = "umi_defined_ids";
    private static final String UMI_HIGHLIGHT = "umi_highlight";
    private static final String UMI_BASE_DIFF_STATS = "umi_base_diff_stats";

    public static final char READ_ID_DELIM = ':';
    public static final String READ_ID_DELIM_STR = String.valueOf(READ_ID_DELIM);

    public UmiConfig(boolean enabled, boolean duplex, final String duplexDelim, boolean highlight, boolean baseStats)
    {
        Enabled = enabled;
        Duplex = duplex;
        DuplexDelim = duplexDelim;
        CollapseReversed = true;
        HighlightConsensus = highlight;
        BaseStats = baseStats;
        PermittedBaseDiff = DEFAULT_MAX_UMI_BASE_DIFF;
        mUmiLength = 0;
        mDefinedUmis = Sets.newHashSet();
    }

    public boolean hasDefinedUmis() { return !mDefinedUmis.isEmpty(); }
    public void addDefinedUmis(final Set<String> umis) { mDefinedUmis.addAll(umis); }

    public static UmiConfig from(final CommandLine cmd)
    {
        UmiConfig umiConfig = new UmiConfig(
                cmd.hasOption(UMI_ENABLED),
                cmd.hasOption(UMI_DUPLEX),
                cmd.getOptionValue(UMI_DUPLEX_DELIM, String.valueOf(Constants.DEFAULT_DUPLEX_UMI_DELIM)),
                cmd.hasOption(UMI_HIGHLIGHT),
                cmd.hasOption(UMI_BASE_DIFF_STATS));

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
            if(!exceedsUmiIdDiff(umiId, definedId, PermittedBaseDiff))
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
        options.addOption(UMI_DUPLEX, false, "UMI duplex enabled");
        options.addOption(UMI_DEFINED_IDS, true, "Optional set of defined UMI IDs in file");
        options.addOption(UMI_HIGHLIGHT, false, "Set consensus read to map-qual 0 to highlight in IGV");
        options.addOption(UMI_BASE_DIFF_STATS, false, "Record base difference stats");
    }
}
