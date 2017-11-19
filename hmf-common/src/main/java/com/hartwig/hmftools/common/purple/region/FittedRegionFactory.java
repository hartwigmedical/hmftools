package com.hartwig.hmftools.common.purple.region;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.BAFUtils;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public class FittedRegionFactory {

    private final Gender gender;
    private final int maxPloidy;
    private final double cnvRatioWeightFactor;
    private final boolean ploidyPenaltyExperiment;
    private final double observedBafExponent;

    public FittedRegionFactory(final Gender gender, final int maxPloidy, final double cnvRatioWeightFactor,
            final boolean ploidyPenaltyExperiment, final double observedBafExponent) {
        this.gender = gender;
        this.maxPloidy = maxPloidy;
        this.cnvRatioWeightFactor = cnvRatioWeightFactor;
        this.ploidyPenaltyExperiment = ploidyPenaltyExperiment;
        this.observedBafExponent = observedBafExponent;
    }

    @NotNull
    public List<FittedRegion> fitRegion(final double purity, final double normFactor,
            @NotNull final Collection<ObservedRegion> observedRegions) {

        final Predicate<ObservedRegion> valid = observedRegion -> gender == Gender.MALE || !observedRegion.chromosome().equals("Y");
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @NotNull
    public FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, normFactor);

        double minDeviation = 0;
        double observedBAF = observedRegion.observedBAF();
        double observedTumorRatio = observedRegion.observedTumorRatio();
        double purityAdjustedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());
        double tumorCopyNumber = tumorCopyNumber(observedRegion, purityAdjustedCopyNumber, refNormalisedCopyNumber);
        double tumorBAF = purityAdjuster.purityAdjustedBAF(observedRegion.chromosome(), tumorCopyNumber, observedBAF);

        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                .from(observedRegion)
                .broadBAF(0)
                .broadTumorCopyNumber(0)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .tumorCopyNumber(tumorCopyNumber)
                .tumorBAF(tumorBAF)
                .refNormalisedCopyNumber(Doubles.replaceNaNWithZero(refNormalisedCopyNumber));

        for (int ploidy = 0; ploidy <= maxPloidy; ploidy++) {
            double modelRatio = modelRatio(purity, normFactor, ploidy);
            double cnvDeviation = cnvDeviation(cnvRatioWeightFactor, modelRatio, observedTumorRatio);

            double[] modelBAFWithDeviation =
                    observedRegion.bafCount() == 0 ? new double[] { 0, 0, 0 } : modelBAFToMinimizeDeviation(purity, ploidy, observedBAF);

            double modelBAF = modelBAFWithDeviation[0];
            double bafDeviation = modelBAFWithDeviation[1];

            double ploidyPenalty = PloidyPenalty.penalty(ploidy, (int) modelBAFWithDeviation[2]);
            double deviation = ploidyPenalty * (bafDeviation + cnvDeviation) * Math.pow(observedBAF, observedBafExponent);

            if (ploidy == 1 || deviation < minDeviation) {
                builder.modelPloidy(ploidy)
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

    private double tumorCopyNumber(@NotNull final ObservedRegion region, double purityAdjustedCopyNumber, double refNormalisedCopyNumber) {
        switch (region.status()) {
            case GERMLINE_NOISE:
            case GERMLINE_HOM_DELETION:
            case GERMLINE_HET_DELETION:
            case GERMLINE_AMPLIFICATION:
                return refNormalisedCopyNumber;
        }

        return purityAdjustedCopyNumber;
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
        double finalBAF = 0;
        double finalDeviation = 0;
        int finalMajorAllele = 0;

        int minBetaAllele = BAFUtils.minAlleleCount(ploidy);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            double modelBAF = BAFUtils.modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = bafDeviation(modelBAF, actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < finalDeviation) {
                finalBAF = modelBAF;
                finalDeviation = modelDeviation;
                finalMajorAllele = betaAllele;
            }
        }

        return new double[] { finalBAF, finalDeviation, finalMajorAllele };
    }
}
