package com.hartwig.hmftools.sage.quality;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class QualityConfig
{
    public final double JitterPenalty;
    public final int JitterMinRepeatCount;
    public final int BaseQualityFixedPenalty;
    public final int DistanceFromReadEdgeFixedPenalty;
    public final int DistanceFromReadEdgeFactor;
    public final int FixedPenalty;
    public final double ReadEventsPenalty;
    public final int ImproperPairPenalty;

    private static final String JITTER_PENALTY = "jitter_penalty";
    private static final String JITTER_MIN_REPEAT_COUNT = "jitter_min_repeat_count";
    private static final String BASE_QUAL_FIXED_PENALTY = "base_qual_fixed_penalty";
    private static final String READ_EDGE_FIXED_PENALTY = "read_edge_fixed_penalty";
    private static final String READ_EDGE_FACTOR = "read_edge_factor";
    private static final String MAP_QUAL_FIXED_PENALTY = "fixed_qual_penalty";
    private static final String MAP_QUAL_IMPROPER_PAIR_PENALTY = "improper_pair_qual_penalty";
    private static final String MAP_QUAL_READ_EVENTS_PENALTY = "read_events_qual_penalty";

    private static final double DEFAULT_JITTER_PENALTY = 0.25;
    private static final int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    private static final int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;
    private static final int DEFAULT_READ_EDGE_FIXED_PENALTY = 0;
    private static final int DEFAULT_MAP_QUAL_FIXED_PENALTY = 15;
    private static final int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    private static final double DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY = 7;

    public QualityConfig(final ConfigBuilder configBuilder)
    {
        JitterPenalty = configBuilder.getDecimal(JITTER_PENALTY);
        JitterMinRepeatCount = configBuilder.getInteger(JITTER_MIN_REPEAT_COUNT);
        BaseQualityFixedPenalty = configBuilder.getInteger(BASE_QUAL_FIXED_PENALTY);
        DistanceFromReadEdgeFixedPenalty = configBuilder.getInteger(READ_EDGE_FIXED_PENALTY);
        DistanceFromReadEdgeFactor = configBuilder.getInteger(READ_EDGE_FACTOR);
        FixedPenalty = configBuilder.getInteger(MAP_QUAL_FIXED_PENALTY);
        ReadEventsPenalty = configBuilder.getDecimal(MAP_QUAL_READ_EVENTS_PENALTY);
        ImproperPairPenalty = configBuilder.getInteger(MAP_QUAL_IMPROPER_PAIR_PENALTY);
    }

    public QualityConfig()
    {
        JitterPenalty = DEFAULT_JITTER_PENALTY;
        JitterMinRepeatCount = DEFAULT_JITTER_MIN_REPEAT_COUNT;
        BaseQualityFixedPenalty = DEFAULT_BASE_QUAL_FIXED_PENALTY;
        DistanceFromReadEdgeFixedPenalty = DEFAULT_READ_EDGE_FIXED_PENALTY;
        DistanceFromReadEdgeFactor = 0;
        FixedPenalty = DEFAULT_MAP_QUAL_FIXED_PENALTY;
        ReadEventsPenalty = DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY;
        ImproperPairPenalty = DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY;
    }

    public boolean isHighlyPolymorphic(final GenomePosition position)
    {
        return HlaCommon.containsPosition(position);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(
                JITTER_PENALTY, "Penalty to apply to qual score when read context matches with jitter", DEFAULT_JITTER_PENALTY);

        configBuilder.addInteger(
                JITTER_MIN_REPEAT_COUNT,"Minimum repeat count before applying jitter penalty", DEFAULT_JITTER_MIN_REPEAT_COUNT);

        configBuilder.addInteger(
                BASE_QUAL_FIXED_PENALTY, "Fixed penalty to apply to base quality", DEFAULT_BASE_QUAL_FIXED_PENALTY);

        configBuilder.addInteger(
                MAP_QUAL_FIXED_PENALTY,  "Fixed penalty to apply to map quality", DEFAULT_MAP_QUAL_FIXED_PENALTY);

        configBuilder.addInteger(
                READ_EDGE_FIXED_PENALTY, "Fixed penalty to apply to distance from read edge", DEFAULT_READ_EDGE_FIXED_PENALTY);

        configBuilder.addInteger(READ_EDGE_FACTOR, "Distance from read edge factor", 0);

        configBuilder.addInteger(
                MAP_QUAL_IMPROPER_PAIR_PENALTY,
                "Penalty to apply to map qual when SAM record does not have the ProperPair flag",
                DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY);

        configBuilder.addDecimal(
                MAP_QUAL_READ_EVENTS_PENALTY,
                "Penalty to apply to map qual for additional distance from ref", DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY);
    }
}
