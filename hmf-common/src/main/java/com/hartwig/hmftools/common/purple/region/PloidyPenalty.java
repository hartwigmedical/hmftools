package com.hartwig.hmftools.common.purple.region;

class PloidyPenalty {

    static double penalty(double ploidyPenaltyFactor, double majorAllele, double minorAllele) {
        double wholeGenomeDoublingDistance = 1 + (Math.abs(majorAllele - 2)) + (Math.abs(minorAllele - 2));
        double singleEventDistance = (Math.abs(majorAllele - 1)) + (Math.abs(minorAllele - 1));

        return 1 + ploidyPenaltyFactor * Math.min(singleEventDistance, wholeGenomeDoublingDistance);
    }

}
