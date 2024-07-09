package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BASE_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HIGH_DEPTH_BASE_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HIGH_DEPTH_MAP_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HIGH_DEPTH_MAP_QUAL_RATIO_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_JITTER_MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAP_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MQ_RATIO_FACTOR;

import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class QualityConfig
{
    public final int JitterMinRepeatCount;
    public final int BaseQualityFixedPenalty;
    public final int FixedMapQualPenalty;
    public final double ReadMapQualEventsPenalty;
    public final int ImproperPairPenalty;
    public final double MapQualityRatioFactor;
    public final boolean HighDepthMode;
    public final int HighBaseQualLimit;

    private static final String JITTER_MIN_REPEAT_COUNT = "jitter_min_repeat_count";
    private static final String BASE_QUAL_FIXED_PENALTY = "base_qual_fixed_penalty";
    private static final String MAP_QUAL_FIXED_PENALTY = "fixed_qual_penalty";
    private static final String MAP_QUAL_IMPROPER_PAIR_PENALTY = "improper_pair_qual_penalty";
    private static final String MAP_QUAL_READ_EVENTS_PENALTY = "read_events_qual_penalty";
    private static final String MAP_QUAL_RATIO_FACTOR = "map_qual_ratio_factor";
    private static final String HIGH_DEPTH_BASE_QUAL_LIMIT = "high_depth_base_qual";
    private static final String HIGH_DEPTH_MODE = "high_depth_mode";

    public QualityConfig(final ConfigBuilder configBuilder)
    {
        JitterMinRepeatCount = configBuilder.getInteger(JITTER_MIN_REPEAT_COUNT);
        BaseQualityFixedPenalty = configBuilder.getInteger(BASE_QUAL_FIXED_PENALTY);
        ReadMapQualEventsPenalty = configBuilder.getDecimal(MAP_QUAL_READ_EVENTS_PENALTY);
        ImproperPairPenalty = configBuilder.getInteger(MAP_QUAL_IMPROPER_PAIR_PENALTY);

        HighDepthMode = configBuilder.hasFlag(HIGH_DEPTH_MODE);

        FixedMapQualPenalty = HighDepthMode && !configBuilder.hasValue(MAP_QUAL_FIXED_PENALTY) ?
                DEFAULT_HIGH_DEPTH_MAP_QUAL_FIXED_PENALTY : configBuilder.getInteger(MAP_QUAL_FIXED_PENALTY);

        MapQualityRatioFactor = HighDepthMode && !configBuilder.hasValue(MAP_QUAL_RATIO_FACTOR) ?
                DEFAULT_HIGH_DEPTH_MAP_QUAL_RATIO_FACTOR : configBuilder.getDecimal(MAP_QUAL_RATIO_FACTOR);

        HighBaseQualLimit = HighDepthMode ? configBuilder.getInteger(HIGH_DEPTH_BASE_QUAL_LIMIT) : 0;
    }

    public QualityConfig(boolean highDepthMode)
    {
        JitterMinRepeatCount = DEFAULT_JITTER_MIN_REPEAT_COUNT;
        BaseQualityFixedPenalty = DEFAULT_BASE_QUAL_FIXED_PENALTY;
        FixedMapQualPenalty = DEFAULT_MAP_QUAL_FIXED_PENALTY;
        ReadMapQualEventsPenalty = DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY;
        ImproperPairPenalty = DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY;
        MapQualityRatioFactor = DEFAULT_MQ_RATIO_FACTOR;

        HighDepthMode = highDepthMode;
        HighBaseQualLimit = HighDepthMode ? DEFAULT_HIGH_DEPTH_BASE_QUAL : 0;
    }

    public boolean isHighlyPolymorphic(final BasePosition position)
    {
        return HlaCommon.containsPosition(position);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(
                JITTER_MIN_REPEAT_COUNT,"Minimum repeat count before applying jitter penalty", DEFAULT_JITTER_MIN_REPEAT_COUNT);

        configBuilder.addInteger(
                BASE_QUAL_FIXED_PENALTY, "Fixed penalty to apply to base quality", DEFAULT_BASE_QUAL_FIXED_PENALTY);

        configBuilder.addInteger(
                MAP_QUAL_FIXED_PENALTY,  "Fixed penalty to apply to map quality", DEFAULT_MAP_QUAL_FIXED_PENALTY);

        configBuilder.addInteger(
                MAP_QUAL_IMPROPER_PAIR_PENALTY,
                "Penalty to apply to map qual when SAM record does not have the ProperPair flag",
                DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY);

        configBuilder.addDecimal(
                MAP_QUAL_READ_EVENTS_PENALTY,
                "Penalty to apply to map qual for additional distance from ref", DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY);

        configBuilder.addDecimal(MAP_QUAL_RATIO_FACTOR, "Map quality ratio factor (0 = disabled)", DEFAULT_MQ_RATIO_FACTOR);

        configBuilder.addFlag(HIGH_DEPTH_MODE, "Enable high-depth mode");
        configBuilder.addInteger(HIGH_DEPTH_BASE_QUAL_LIMIT, "High-depth mode min base qual", DEFAULT_HIGH_DEPTH_BASE_QUAL);
    }
}
