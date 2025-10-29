package com.hartwig.hmftools.redux.bqr;

import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BqrConfig
{
    public final boolean Enabled;
    public final boolean UseAllRegions;
    public final boolean WritePlot;

    private static final String DISABLE_BQR = "bqr_disable";
    private static final String WRITE_BQR_PLOT = "bqr_write_plot";
    private static final String USE_ALL_REGIONS = "bqr_use_all_regions";

    private static final String CFG_BQR_LOG_TNCS = "bqr_debug_tnc";
    private static final String CFG_BQR_LOG_QUALS = "bqr_debug_qual";
    private static final String CFG_BQR_LOG_POSITIONS = "bqr_debug_position";
    private static final String CFG_BQR_LOG_CON_TYPES = "bqr_debug_consensus";

    // debug only
    protected static boolean LogDebug = false;
    protected static final List<String> LOG_TNCS = Lists.newArrayList(); // eg AAT;AGT
    protected static final List<Byte> LOG_QUAL = Lists.newArrayList(); // eg 25;35
    protected static final List<Integer> LOG_POSITIONS = Lists.newArrayList();
    protected static final List<ConsensusType> LOG_CONSENSUS_TYPES = Lists.newArrayList();

    public BqrConfig(final ConfigBuilder configBuilder)
    {
        Enabled = !configBuilder.hasFlag(DISABLE_BQR);
        WritePlot = configBuilder.hasFlag(WRITE_BQR_PLOT);
        UseAllRegions = configBuilder.hasFlag(USE_ALL_REGIONS);

        if(configBuilder.hasValue(CFG_BQR_LOG_TNCS))
        {
            Arrays.stream(configBuilder.getValue(CFG_BQR_LOG_TNCS).split(ITEM_DELIM, -1)).forEach(x -> LOG_TNCS.add(x));
        }

        if(configBuilder.hasValue(CFG_BQR_LOG_POSITIONS))
        {
            Arrays.stream(configBuilder.getValue(CFG_BQR_LOG_POSITIONS).split(ITEM_DELIM, -1))
                    .forEach(x -> LOG_POSITIONS.add(Integer.parseInt(x)));
        }

        if(configBuilder.hasValue(CFG_BQR_LOG_QUALS))
        {
            Arrays.stream(configBuilder.getValue(CFG_BQR_LOG_QUALS).split(ITEM_DELIM, -1)).forEach(x -> LOG_QUAL.add(Byte.parseByte(x)));
        }

        if(configBuilder.hasValue(CFG_BQR_LOG_CON_TYPES))
        {
            Arrays.stream(configBuilder.getValue(CFG_BQR_LOG_CON_TYPES).split(ITEM_DELIM, -1))
                    .forEach(x -> LOG_CONSENSUS_TYPES.add(ConsensusType.valueOf(x)));
        }

        LogDebug = !LOG_TNCS.isEmpty() || !LOG_POSITIONS.isEmpty() || !LOG_QUAL.isEmpty() || !LOG_CONSENSUS_TYPES.isEmpty();
    }

    public BqrConfig()
    {
        Enabled = false;
        WritePlot = false;
        UseAllRegions = false;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(DISABLE_BQR, "Disable Base Quality Recalibration");
        configBuilder.addFlag(WRITE_BQR_PLOT, "Generate BQR plot");
        configBuilder.addFlag(USE_ALL_REGIONS, "Run BQR from full BAM");

        configBuilder.addConfigItem(CFG_BQR_LOG_TNCS, "BQR log debug tri-nuc contexts, separated by ';'");
        configBuilder.addConfigItem(CFG_BQR_LOG_QUALS, "BQR log debug quals, separated by ';'");
        configBuilder.addConfigItem(CFG_BQR_LOG_POSITIONS, "BQR log debug positions, separated by ';'");
        configBuilder.addConfigItem(CFG_BQR_LOG_CON_TYPES, "BQR log debug consensus types, separated by ';'");
    }
}
