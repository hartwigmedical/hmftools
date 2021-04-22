package com.hartwig.hmftools.purple.region;

final class EventPenalty {

    static double penalty(double eventPenaltyFactor, double majorAllele, double minorAllele) {
        double wholeGenomeDoublingDistance = 1 + (Math.abs(majorAllele - 2)) + (Math.abs(minorAllele - 2));
        double singleEventDistance = (Math.abs(majorAllele - 1)) + (Math.abs(minorAllele - 1));

        return 1 + eventPenaltyFactor * Math.min(singleEventDistance, wholeGenomeDoublingDistance);
    }
}
