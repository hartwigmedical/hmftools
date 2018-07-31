package com.hartwig.hmftools.common.purple.region;

class PloidyPenalty {

    static double penalty(int ploidy, int majorAllele) {
        int minorAllele = ploidy - majorAllele;

        int wholeGenomeDoublingDistance = 1 + Math.abs(majorAllele - 2) + Math.abs(minorAllele - 2);
        int singleEventDistance = Math.abs(majorAllele - 1) + Math.abs(minorAllele - 1);

        return 1 + Math.min(singleEventDistance, wholeGenomeDoublingDistance);
    }

    static double penaltyv2(double ploidyPenaltyFactor, double majorAllele, double minorAllele) {
        double wholeGenomeDoublingDistance = 1 + (Math.abs(majorAllele - 2)) + (Math.abs(minorAllele - 2));
        double singleEventDistance = (Math.abs(majorAllele - 1)) + (Math.abs(minorAllele - 1));

        return 1 + ploidyPenaltyFactor * Math.min(singleEventDistance, wholeGenomeDoublingDistance);
    }

}
