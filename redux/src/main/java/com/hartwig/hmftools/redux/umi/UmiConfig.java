package com.hartwig.hmftools.redux.umi;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_MAX_UMI_BASE_DIFF;
import static com.hartwig.hmftools.redux.umi.UmiUtils.exceedsUmiIdDiff;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.common.Constants;

public class UmiConfig
{
    public final boolean Enabled;
    public final boolean Duplex; // collapse duplex UMI groups
    public final String DuplexDelim;
    public final boolean CollapseReversed;
    public final boolean BaseStats;
    public final int PermittedBaseDiff;

    private int mUmiLength; // set and accessed in a thread-safe way

    private final Set<String> mDefinedUmis;

    // config options
    private static final String UMI_ENABLED = "umi_enabled";
    private static final String UMI_DUPLEX = "umi_duplex";
    private static final String UMI_DUPLEX_DELIM = "umi_duplex_delim";
    private static final String UMI_DEFINED_IDS = "umi_defined_ids";
    private static final String UMI_BASE_DIFF_STATS = "umi_base_diff_stats";

    public static final char READ_ID_DELIM = ':';
    public static final String READ_ID_DELIM_STR = String.valueOf(READ_ID_DELIM);

    public UmiConfig(boolean enabled, boolean duplex, final String duplexDelim, boolean baseStats)
    {
        Enabled = enabled;
        Duplex = duplex;
        DuplexDelim = duplexDelim;
        CollapseReversed = true;
        BaseStats = baseStats;
        PermittedBaseDiff = DEFAULT_MAX_UMI_BASE_DIFF;
        mUmiLength = 0;
        mDefinedUmis = Sets.newHashSet();
    }

    public boolean hasDefinedUmis() { return !mDefinedUmis.isEmpty(); }
    public void addDefinedUmis(final Set<String> umis) { mDefinedUmis.addAll(umis); }

    public static UmiConfig from(final ConfigBuilder configBuilder)
    {
        UmiConfig umiConfig = new UmiConfig(
                configBuilder.hasFlag(UMI_ENABLED),
                configBuilder.hasFlag(UMI_DUPLEX),
                configBuilder.getValue(UMI_DUPLEX_DELIM),
                configBuilder.hasFlag(UMI_BASE_DIFF_STATS));

        if(configBuilder.hasValue(UMI_DEFINED_IDS))
        {
            String definedUmiIdsFilename = configBuilder.getValue(UMI_DEFINED_IDS);
            try
            {
                List<String> umis = Files.readAllLines(Paths.get(definedUmiIdsFilename));
                RD_LOGGER.info("loaded {} defined UMI IDs from {}", umis.size(), definedUmiIdsFilename);
                umiConfig.addDefinedUmis(umis.stream().collect(Collectors.toSet()));
            }
            catch(IOException e)
            {
                RD_LOGGER.error("failed to load defined UMI IDs: {}", e.toString());
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

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(UMI_ENABLED, "Use UMIs for duplicates");
        configBuilder.addFlag(UMI_DUPLEX, "UMI duplex enabled");
        configBuilder.addPath(UMI_DEFINED_IDS, false, "Optional set of defined UMI IDs in file");
        configBuilder.addFlag(UMI_BASE_DIFF_STATS, "Record base difference stats");

        configBuilder.addConfigItem(
                UMI_DUPLEX_DELIM, false,
                "UMI duplex delimiter, default: " + Constants.DEFAULT_DUPLEX_UMI_DELIM,
                String.valueOf(Constants.DEFAULT_DUPLEX_UMI_DELIM));
    }
}
