package com.hartwig.hmftools.common.variant;

public enum Clonality {
    CLONAL,
    SUBCLONAL,
    INCONSISTENT,
    UNKNOWN;

    private static final int MULTIPLIER = 10_000_000;

//    public static Clonality fromSample(double purity, double normFactor, double observedTumorRatio, double obsevedNormalRatio,
//            @NotNull AllelicDepth depth) {
//
//        double monoploidProbability = purity * normFactor / 2 / observedTumorRatio;
//        final BinomialDistribution monoploidDistribution =
//                new BinomialDistribution(depth.totalReadCount() * MULTIPLIER, monoploidProbability / MULTIPLIER);
//        if (depth.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.01)) {
//            return SUBCLONAL;
//        }
//
//        double inconsistentProbability = (purity - 1) * obsevedNormalRatio * normFactor / observedTumorRatio + 1;
//        final BinomialDistribution inconsistentDistribution =
//                new BinomialDistribution(depth.totalReadCount() * MULTIPLIER, inconsistentProbability / MULTIPLIER);
//        if (depth.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
//            return INCONSISTENT;
//        }
//
//        return CLONAL;
//    }
}
