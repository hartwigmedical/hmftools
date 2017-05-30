package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedBAF;
import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedCopyNumber;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public class FittedRegionFactory {

    @VisibleForTesting
    static final double NORMAL_BAF = 0.533;

    private final int maxPloidy;
    private final double cnvRatioWeightFactor;

    public FittedRegionFactory(final int maxPloidy, final double cnvRatioWeightFactor) {
        this.maxPloidy = maxPloidy;
        this.cnvRatioWeightFactor = cnvRatioWeightFactor;
    }

    @NotNull
    public List<FittedRegion> fitRegion(final double purity, final double normFactor,
            @NotNull final Collection<ObservedRegion> observedRegions) {
        return observedRegions.stream().map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @NotNull
    public FittedRegion fitRegion(final double purity, final double normFactor,
            final @NotNull ObservedRegion observedRegion) {
        double minDeviation = 0;
        double observedBAF = observedRegion.observedBAF();
        double observedTumorRatio = observedRegion.observedTumorRatio();
        double tumorCopyNumber = purityAdjustedCopyNumber(purity, normFactor, observedTumorRatio);

        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder().from(observedRegion).status(
                FreecStatus.fromNormalRatio(observedRegion.observedNormalRatio())).broadBAF(0).broadTumorCopyNumber(
                0).segmentBAF(0).segmentTumorCopyNumber(0).tumorCopyNumber(tumorCopyNumber).refNormalisedCopyNumber(
                Doubles.replaceNaNWithZero(
                        observedTumorRatio / observedRegion.observedNormalRatio() / normFactor * 2));

        for (int ploidy = 1; ploidy <= maxPloidy; ploidy++) {
            double modelRatio = modelRatio(purity, normFactor, ploidy);
            double cnvDeviation = cnvDeviation(cnvRatioWeightFactor, modelRatio, observedTumorRatio);

            double[] modelBAFWithDeviation = observedRegion.bafCount() == 0 ?
                    new double[] { 0, 0 } :
                    modelBAFToMinimizeDeviation(purity, ploidy, observedBAF);

            double modelBAF = modelBAFWithDeviation[0];
            double bafDeviation = modelBAFWithDeviation[1];

            double deviation =
                    Math.pow(Math.max(ploidy, 1.5) / 2.0, 0.85) * (bafDeviation + cnvDeviation) * observedBAF;

            if (ploidy == 1 || deviation < minDeviation) {
                builder.fittedPloidy(ploidy).modelBAF(modelBAF).modelTumorRatio(modelRatio).bafDeviation(
                        bafDeviation).cnvDeviation(cnvDeviation).purityAdjustedBAF(
                        purityAdjustedBAF(purity, ploidy, observedBAF)).deviation(deviation);
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
    static double cnvDeviation(final double cnvRatioWeighFactor, final double modelCNVRatio,
            final double actualRatio) {
        return cnvRatioWeighFactor * Math.abs(modelCNVRatio - actualRatio);
    }

    @VisibleForTesting
    static double bafDeviation(final boolean heterozygous, final double modelBAF, final double actualBAF) {
        if (heterozygous && Doubles.lessOrEqual(actualBAF, NORMAL_BAF)) {
            return 0;
        }

        return Math.abs(modelBAF - actualBAF);
    }

    @VisibleForTesting
    static double[] modelBAFToMinimizeDeviation(final double purity, final int ploidy, final double actualBAF) {
        double result = 0;
        double deviation = 0;

        int minBetaAllele = (int) Math.round(ploidy / 2d);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            boolean isHeterozygous = ploidy / betaAllele == 2;
            double modelBAF = isHeterozygous ? NORMAL_BAF : modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = bafDeviation(isHeterozygous, modelBAF, actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < deviation) {
                result = modelBAF;
                deviation = modelDeviation;
            }
        }

        return new double[] { result, deviation };
    }

    @VisibleForTesting
    public static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        if (ploidy / alleleCount == 2) {
            return NORMAL_BAF;
        }

        return (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2));
    }
}
