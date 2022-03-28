package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.hla.HlaCommon;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class QualityConfig
{
    public final double JitterPenalty;
    public final int JitterMinRepeatCount;
    public final int BaseQualityFixedPenalty;
    public final int DistanceFromReadEdgeFixedPenalty;
    public final int MapQualityFixedPenalty;
    public final int MapQualityReadEventsPenalty;
    public final int MapQualityImproperPairPenalty;

    private static final String JITTER_PENALTY = "jitter_penalty";
    private static final String JITTER_MIN_REPEAT_COUNT = "jitter_min_repeat_count";
    private static final String BASE_QUAL_FIXED_PENALTY = "base_qual_fixed_penalty";
    private static final String READ_EDGE_FIXED_PENALTY = "read_edge_fixed_penalty";
    private static final String MAP_QUAL_FIXED_PENALTY = "map_qual_fixed_penalty";
    private static final String MAP_QUAL_IMPROPER_PAIR_PENALTY = "map_qual_improper_pair_penalty";
    private static final String MAP_QUAL_READ_EVENTS_PENALTY = "map_qual_read_events_penalty";

    private static final double DEFAULT_JITTER_PENALTY = 0.25;
    private static final int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    private static final int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;
    private static final int DEFAULT_READ_EDGE_FIXED_PENALTY = 0;
    private static final int DEFAULT_MAP_QUAL_FIXED_PENALTY = 15;
    private static final int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    private static final int DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY = 8;

    public QualityConfig(final CommandLine cmd)
    {
        JitterPenalty = getConfigValue(cmd, JITTER_PENALTY, DEFAULT_JITTER_PENALTY);
        JitterMinRepeatCount = getConfigValue(cmd, JITTER_MIN_REPEAT_COUNT, DEFAULT_JITTER_MIN_REPEAT_COUNT);
        BaseQualityFixedPenalty = getConfigValue(cmd, BASE_QUAL_FIXED_PENALTY, DEFAULT_BASE_QUAL_FIXED_PENALTY);
        DistanceFromReadEdgeFixedPenalty = getConfigValue(cmd, READ_EDGE_FIXED_PENALTY, DEFAULT_READ_EDGE_FIXED_PENALTY);
        MapQualityFixedPenalty = getConfigValue(cmd, MAP_QUAL_FIXED_PENALTY, DEFAULT_MAP_QUAL_FIXED_PENALTY);
        MapQualityReadEventsPenalty = getConfigValue(cmd, MAP_QUAL_READ_EVENTS_PENALTY, DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY);
        MapQualityImproperPairPenalty = getConfigValue(cmd, MAP_QUAL_IMPROPER_PAIR_PENALTY, DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY);
    }

    public QualityConfig()
    {
        JitterPenalty = DEFAULT_JITTER_PENALTY;
        JitterMinRepeatCount = DEFAULT_JITTER_MIN_REPEAT_COUNT;
        BaseQualityFixedPenalty = DEFAULT_BASE_QUAL_FIXED_PENALTY;
        DistanceFromReadEdgeFixedPenalty = DEFAULT_READ_EDGE_FIXED_PENALTY;
        MapQualityFixedPenalty = DEFAULT_MAP_QUAL_FIXED_PENALTY;
        MapQualityReadEventsPenalty = DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY;
        MapQualityImproperPairPenalty = DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY;
    }

    public boolean isHighlyPolymorphic(final GenomePosition position)
    {
        return HlaCommon.containsPosition(position);
    }

    public static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(JITTER_PENALTY,
                true,
                "Penalty to apply to qual score when read context matches with jitter [" + DEFAULT_JITTER_PENALTY + "]");
        options.addOption(JITTER_MIN_REPEAT_COUNT,
                true,
                "Minimum repeat count before applying jitter penalty [" + DEFAULT_JITTER_MIN_REPEAT_COUNT + "]");
        options.addOption(BASE_QUAL_FIXED_PENALTY,
                true,
                "Fixed penalty to apply to base quality [" + DEFAULT_BASE_QUAL_FIXED_PENALTY + "]");
        options.addOption(MAP_QUAL_FIXED_PENALTY, true, "Fixed penalty to apply to map quality [" + DEFAULT_MAP_QUAL_FIXED_PENALTY + "]");
        options.addOption(READ_EDGE_FIXED_PENALTY,
                true,
                "Fixed penalty to apply to distance from read edge [" + DEFAULT_READ_EDGE_FIXED_PENALTY + "]");
        options.addOption(MAP_QUAL_IMPROPER_PAIR_PENALTY,
                true,
                "Penalty to apply to map qual when SAM record does not have the ProperPair flag [" + DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY
                        + "]");
        options.addOption(MAP_QUAL_READ_EVENTS_PENALTY,
                true,
                "Penalty to apply to map qual for additional distance from ref [" + DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY + "]");

        return options;
    }

}
