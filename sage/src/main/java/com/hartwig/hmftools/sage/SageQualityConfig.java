package com.hartwig.hmftools.sage;

public interface SageQualityConfig {

    default int baseQualityFixedPenalty() {
        return 12;
    }

    default int mapQualityFixedPenalty() {
        return 4;
    }

    default int mapQualityAdditionalDistanceFromRefPenalty() {
        return 4;
    }

    default int mapQualityImproperPairPenalty() {
        return 5;
    }
}
