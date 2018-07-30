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

public class FittedRegionFactoryV1 implements FittedRegionFactory {

    private final Gender gender;
    private final int maxPloidy;
    private final double cnvRatioWeightFactor;
    private final double observedBafExponent;
    private final BAFUtils bafUtils;

    public FittedRegionFactoryV1(final Gender gender, final int maxPloidy, final double cnvRatioWeightFactor,
            final int averageReadDepth, final double observedBafExponent) {
        this.gender = gender;
        this.maxPloidy = maxPloidy;
        this.cnvRatioWeightFactor = cnvRatioWeightFactor;
        this.observedBafExponent = observedBafExponent;
        bafUtils = new BAFUtils(averageReadDepth);
    }

    @Override
    @NotNull
    public List<FittedRegion> fitRegion(final double purity, final double normFactor,
            @NotNull final Collection<ObservedRegion> observedRegions) {

        final Predicate<ObservedRegion> valid = observedRegion -> !(gender == Gender.FEMALE && observedRegion.chromosome().equals("Y"));
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @Override
    @NotNull
    public FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, normFactor);

        double minDeviation = 0;
        double observedBAF = observedRegion.observedBAF();
        double observedTumorRatio = observedRegion.observedTumorRatio();
        double purityAdjustedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());
        double tumorBAF = bafUtils.purityAdjustedBAF(purityAdjuster, observedRegion.chromosome(), purityAdjustedCopyNumber, observedBAF);

        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                .from(observedRegion)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .tumorCopyNumber(purityAdjustedCopyNumber)
                .tumorBAF(tumorBAF)
                .refNormalisedCopyNumber(Doubles.replaceNaNWithZero(refNormalisedCopyNumber));

        for (int ploidy = 0; ploidy <= maxPloidy; ploidy++) {
            double modelRatio = modelRatio(purity, normFactor, ploidy);
            double cnvDeviation = cnvDeviation(cnvRatioWeightFactor, modelRatio, observedTumorRatio);

            double[] modelBAFWithDeviation =
                    observedRegion.bafCount() == 0 ? new double[] { 0, 0, 0 } : bafUtils.modelBAFToMinimizeDeviation(purity, ploidy, observedBAF);

            double modelBAF = modelBAFWithDeviation[0];
            double bafDeviation = modelBAFWithDeviation[1];

            double ploidyPenalty = PloidyPenalty.penalty(ploidy, (int) modelBAFWithDeviation[2]);
            double deviation = (bafDeviation + cnvDeviation) * Math.pow(observedBAF, observedBafExponent);

            if (ploidy == 1 || deviation < minDeviation) {
                builder.modelPloidy(ploidy)
                        .modelBAF(modelBAF)
                        .modelTumorRatio(modelRatio)
                        .bafDeviation(bafDeviation)
                        .cnvDeviation(cnvDeviation)
                        .ploidyPenalty(ploidyPenalty)
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

}
