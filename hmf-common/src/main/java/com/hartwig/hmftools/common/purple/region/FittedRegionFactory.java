package com.hartwig.hmftools.common.purple.region;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public class FittedRegionFactory {

    public static final double NORMAL_BAF = 0.535;

    private final Gender gender;
    private final int maxPloidy;
    private final double cnvRatioWeightFactor;
    private final double ploidyPenaltyExponent;
    private final double observedBafExponent;

    public FittedRegionFactory(final Gender gender, final int maxPloidy, final double cnvRatioWeightFactor,
            final double ploidyPenaltyExponent, final double observedBafExponent) {
        this.gender = gender;
        this.maxPloidy = maxPloidy;
        this.cnvRatioWeightFactor = cnvRatioWeightFactor;
        this.ploidyPenaltyExponent = ploidyPenaltyExponent;
        this.observedBafExponent = observedBafExponent;
    }

    @NotNull
    public List<FittedRegion> fitRegion(final double purity, final double normFactor, @NotNull final Collection<ObservedRegion> observedRegions) {
        return observedRegions.stream().map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @NotNull
    public FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, normFactor);

        double minDeviation = 0;
        double observedBAF = observedRegion.observedBAF();
        double observedTumorRatio = observedRegion.observedTumorRatio();
        double tumorCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);

        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                .from(observedRegion)
                .status(FreecStatus.fromNormalRatio(observedRegion.observedNormalRatio()))
                .broadBAF(0)
                .broadTumorCopyNumber(0)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(Doubles.replaceNaNWithZero(observedTumorRatio / observedRegion.observedNormalRatio() / normFactor * 2));

        for (int ploidy = 1; ploidy <= maxPloidy; ploidy++) {
            double modelRatio = modelRatio(purity, normFactor, ploidy);
            double cnvDeviation = cnvDeviation(cnvRatioWeightFactor, modelRatio, observedTumorRatio);

            double[] modelBAFWithDeviation = observedRegion.bafCount() == 0 ? new double[] { 0, 0 } : modelBAFToMinimizeDeviation(purity, ploidy, observedBAF);

            double modelBAF = modelBAFWithDeviation[0];
            double bafDeviation = modelBAFWithDeviation[1];

            double deviation = Math.pow(Math.max(ploidy, 2) / 2.0, ploidyPenaltyExponent) * (bafDeviation + cnvDeviation) * Math.pow(observedBAF, observedBafExponent);

            if (ploidy == 1 || deviation < minDeviation) {
                builder.fittedPloidy(ploidy)
                        .modelBAF(modelBAF)
                        .modelTumorRatio(modelRatio)
                        .bafDeviation(bafDeviation)
                        .cnvDeviation(cnvDeviation)
                        .deviation(deviation);
                minDeviation = deviation;
            }
        }

        return builder.build();
    }

    @VisibleForTesting
    static double modelRatio(final double purity, final double normFactor, final int ploidy) {
        return normFactor + (ploidy - 2) * purity * normFactor / 2d;
    }

    @VisibleForTesting
    static double cnvDeviation(final double cnvRatioWeighFactor, final double modelCNVRatio, final double actualRatio) {
        return cnvRatioWeighFactor * Math.abs(modelCNVRatio - actualRatio);
    }

    @VisibleForTesting
    static double bafDeviation(final double modelBAF, final double actualBAF) {
        return Math.abs(modelBAF - actualBAF);
    }

    @VisibleForTesting
    static double[] modelBAFToMinimizeDeviation(final double purity, final int ploidy, final double actualBAF) {
        double result = 0;
        double deviation = 0;

        int minBetaAllele = (int) Math.round(ploidy / 2d);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            double modelBAF = modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = bafDeviation(modelBAF, actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < deviation) {
                result = modelBAF;
                deviation = modelDeviation;
            }
        }

        return new double[] { result, deviation };
    }

    @VisibleForTesting
    static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        assert (alleleCount >= ploidy / 2d);

        if (ploidy / alleleCount == 2) {
            return NORMAL_BAF;
        }

        return Math.max(NORMAL_BAF, (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2)));
    }
}
