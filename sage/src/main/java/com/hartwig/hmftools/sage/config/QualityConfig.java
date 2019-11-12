package com.hartwig.hmftools.sage.config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface QualityConfig {

    String JITTER_PENALTY = "jitter_penalty";
    String JITTER_MIN_REPEAT_COUNT = "jitter_min_repeat_count";
    String BASE_QUAL_FIXED_PENALTY = "base_qual_fixed_penalty";
    String MAP_QUAL_FIXED_PENALTY = "map_qual_fixed_penalty";
    String MAP_QUAL_IMPROPER_PAIR_PENALTY = "map_qual_improper_pair_penalty";
    String MAP_QUAL_DISTANCE_FROM_REF = "map_qual_distance_from_ref_penalty";

    double DEFAULT_JITTER_PENALTY = 0.25;
    int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;
    int DEFAULT_MAP_QUAL_FIXED_PENALTY = 24;
    int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    int DEFAULT_MAP_QUAL_DISTANCE_FROM_REF = 5;

    double jitterPenalty();

    int jitterMinRepeatCount();

    int baseQualityFixedPenalty();

    int mapQualityFixedPenalty();

    int mapQualityAdditionalDistanceFromRefPenalty();

    int mapQualityImproperPairPenalty();

    default int modifiedMapQuality(int mapQuality, int distance, boolean properPairFlag) {
        int improperPairPenalty = mapQualityImproperPairPenalty() * (properPairFlag ? 0 : 1);
        int distancePenalty = (distance - 1) * mapQualityAdditionalDistanceFromRefPenalty();
        return mapQuality - mapQualityFixedPenalty() - improperPairPenalty - distancePenalty;
    }

    default int modifiedBaseQuality(int baseQuality, int distanceFromReadEdge) {
        return Math.min(baseQuality - baseQualityFixedPenalty(), distanceFromReadEdge) ;
    }

    default double jitterPenalty(int repeatCount) {
        return (jitterPenalty() * Math.max(0, repeatCount - jitterMinRepeatCount()));
    }

    @NotNull
    static Options createOptions() {
        final Options options = new Options();

        options.addOption(JITTER_PENALTY, true, "Penalty to apply to qual score when read context matches with jitter [" + DEFAULT_JITTER_PENALTY + "]");
        options.addOption(JITTER_MIN_REPEAT_COUNT, true, "Minimum repeat count before applying jitter penalty [" + DEFAULT_JITTER_MIN_REPEAT_COUNT + "]");
        options.addOption(BASE_QUAL_FIXED_PENALTY, true, "Fixed penalty to apply to base quality [" + DEFAULT_BASE_QUAL_FIXED_PENALTY + "]");
        options.addOption(MAP_QUAL_FIXED_PENALTY, true, "Fixed penalty to apply to map quality [" + DEFAULT_MAP_QUAL_FIXED_PENALTY + "]");
        options.addOption(MAP_QUAL_IMPROPER_PAIR_PENALTY, true, "Penalty to apply to map qual when SAM record does not have the ProperPair flag [" + DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY + "]");
        options.addOption(MAP_QUAL_DISTANCE_FROM_REF, true, "Penalty to apply to map qual for additional distance from ref [" + DEFAULT_MAP_QUAL_DISTANCE_FROM_REF + "]");

        return options;
    }

    @NotNull
    static QualityConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {

        return ImmutableQualityConfig.builder()
                .jitterPenalty(SageConfig.defaultValue(cmd, JITTER_PENALTY, DEFAULT_JITTER_PENALTY))
                .jitterMinRepeatCount(SageConfig.defaultIntValue(cmd, JITTER_MIN_REPEAT_COUNT, DEFAULT_JITTER_MIN_REPEAT_COUNT))
                .baseQualityFixedPenalty(SageConfig.defaultIntValue(cmd, BASE_QUAL_FIXED_PENALTY, DEFAULT_BASE_QUAL_FIXED_PENALTY))
                .mapQualityFixedPenalty(SageConfig.defaultIntValue(cmd, MAP_QUAL_FIXED_PENALTY, DEFAULT_MAP_QUAL_FIXED_PENALTY))
                .mapQualityAdditionalDistanceFromRefPenalty(SageConfig.defaultIntValue(cmd, MAP_QUAL_DISTANCE_FROM_REF, DEFAULT_MAP_QUAL_DISTANCE_FROM_REF))
                .mapQualityImproperPairPenalty(SageConfig.defaultIntValue(cmd, MAP_QUAL_IMPROPER_PAIR_PENALTY, DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY))
                .build();

    }

}
