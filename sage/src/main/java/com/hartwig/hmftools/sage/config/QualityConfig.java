package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.cli.Configs.defaultDoubleValue;
import static com.hartwig.hmftools.common.cli.Configs.defaultIntValue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface QualityConfig {

    String JITTER_PENALTY = "jitter_penalty";
    String JITTER_MIN_REPEAT_COUNT = "jitter_min_repeat_count";
    String BASE_QUAL_FIXED_PENALTY = "base_qual_fixed_penalty";
    String READ_EDGE_FIXED_PENALTY = "read_edge_fixed_penalty";
    String MAP_QUAL_FIXED_PENALTY = "map_qual_fixed_penalty";
    String MAP_QUAL_IMPROPER_PAIR_PENALTY = "map_qual_improper_pair_penalty";
    String MAP_QUAL_DISTANCE_FROM_REF = "map_qual_distance_from_ref_penalty";

    double DEFAULT_JITTER_PENALTY = 0.25;
    int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;
    int DEFAULT_READ_EDGE_FIXED_PENALTY = 0;
    int DEFAULT_MAP_QUAL_FIXED_PENALTY = 15;
    int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    int DEFAULT_MAP_QUAL_DISTANCE_FROM_REF = 10;

    double jitterPenalty();

    int jitterMinRepeatCount();

    int baseQualityFixedPenalty();

    int distanceFromReadEdgeFixedPenalty();

    int mapQualityFixedPenalty();

    int mapQualityAdditionalDistanceFromRefPenalty();

    int mapQualityImproperPairPenalty();

    default int modifiedMapQuality(int mapQuality, int distance, boolean properPairFlag) {
        int improperPairPenalty = mapQualityImproperPairPenalty() * (properPairFlag ? 0 : 1);
        int distancePenalty = (distance - 1) * mapQualityAdditionalDistanceFromRefPenalty();
        return mapQuality - mapQualityFixedPenalty() - improperPairPenalty - distancePenalty;
    }

    default double modifiedBaseQuality(double baseQuality, int distanceFromReadEdge) {
        return Math.min(baseQuality - baseQualityFixedPenalty(), 3 * distanceFromReadEdge - distanceFromReadEdgeFixedPenalty());
    }

    default double jitterPenalty(int repeatCount) {
        return (jitterPenalty() * Math.max(0, repeatCount - jitterMinRepeatCount()));
    }

    default int baseQualityRecalibrationMaxAltCount() {
        return 3;
    }

    @NotNull
    static Options createOptions() {
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
        options.addOption(MAP_QUAL_DISTANCE_FROM_REF,
                true,
                "Penalty to apply to map qual for additional distance from ref [" + DEFAULT_MAP_QUAL_DISTANCE_FROM_REF + "]");

        return options;
    }

    @NotNull
    static QualityConfig createConfig(@NotNull final CommandLine cmd) {

        return ImmutableQualityConfig.builder()
                .jitterPenalty(defaultDoubleValue(cmd, JITTER_PENALTY, DEFAULT_JITTER_PENALTY))
                .jitterMinRepeatCount(defaultIntValue(cmd, JITTER_MIN_REPEAT_COUNT, DEFAULT_JITTER_MIN_REPEAT_COUNT))
                .baseQualityFixedPenalty(defaultIntValue(cmd, BASE_QUAL_FIXED_PENALTY, DEFAULT_BASE_QUAL_FIXED_PENALTY))
                .distanceFromReadEdgeFixedPenalty(defaultIntValue(cmd, READ_EDGE_FIXED_PENALTY, DEFAULT_READ_EDGE_FIXED_PENALTY))
                .mapQualityFixedPenalty(defaultIntValue(cmd, MAP_QUAL_FIXED_PENALTY, DEFAULT_MAP_QUAL_FIXED_PENALTY))
                .mapQualityAdditionalDistanceFromRefPenalty(defaultIntValue(cmd,
                        MAP_QUAL_DISTANCE_FROM_REF,
                        DEFAULT_MAP_QUAL_DISTANCE_FROM_REF))
                .mapQualityImproperPairPenalty(defaultIntValue(cmd, MAP_QUAL_IMPROPER_PAIR_PENALTY, DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY))
                .build();

    }

}
