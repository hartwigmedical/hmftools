package com.hartwig.hmftools.common.purple.region;

class PloidyPenalty {

    static double penalty(int ploidy) {
        return Math.max(ploidy, 2) / 2.0;
    }

    static double penalty(int ploidy, int majorAllele) {
        int minorAllele = ploidy - majorAllele;

        int wholeGenomeDoublingDistance = 1 + Math.abs(majorAllele - 2) + Math.abs(minorAllele - 2);
        int singleEventDistance = Math.abs(majorAllele - 1) + Math.abs(minorAllele - 1);

        return 1 + Math.min(singleEventDistance, wholeGenomeDoublingDistance);
    }
}
