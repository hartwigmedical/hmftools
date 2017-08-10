package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;

final class CopyNumberBuilder {

    private static final double MIN_COPY_NUMBER_TOLERANCE = 0.3;
    @VisibleForTesting
    static final double LC_MAX_COPY_NUMBER_TOLERANCE = 1.3;
    @VisibleForTesting
    static final double HC_MAX_COPY_NUMBER_TOLERANCE = 0.7;

    @NotNull
    private final PurityAdjuster purityAdjuster;
    private final CombinedFittedRegion combinedRegion;

    CopyNumberBuilder(boolean highConfidence, @NotNull final PurityAdjuster purityAdjuster, final FittedRegion combinedRegion) {
        this.purityAdjuster = purityAdjuster;
        this.combinedRegion = new CombinedFittedRegion(highConfidence, combinedRegion);
    }

    public String chromosome() {
        return combinedRegion.region().chromosome();
    }

    public final int bafCount() {
        return combinedRegion.region().bafCount();
    }

    public final double averageObservedBAF() {
        return combinedRegion.region().observedBAF();
    }

    public final double averageTumorCopyNumber() {
        return combinedRegion.region().tumorCopyNumber();
    }

    public final double averageRefNormalisedCopyNumber() {
        return combinedRegion.region().refNormalisedCopyNumber();
    }

    public void extendRegion(@NotNull final FittedRegion value) {
        assert (chromosome().equals(value.chromosome())) : "Regions cannot be extended between chromosomes";
        combinedRegion.combine(value);
    }

    @NotNull
    public FittedRegion build() {
        return combinedRegion.region();
    }


    boolean withinCopyNumberTolerance(@NotNull final FittedRegion copyNumber) {
        double tumorCopyNumberDeviation = Math.abs(copyNumber.tumorCopyNumber() - averageTumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(copyNumber.refNormalisedCopyNumber() - averageRefNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);
        double rawMaxDeviation = maxCopyNumberDeviation(copyNumber);
        double adjustedMaxDeviation = purityAdjuster.purityAdjustedMaxCopyNumberDeviation(rawMaxDeviation);

        return Doubles.lessOrEqual(copyNumberDeviation, adjustedMaxDeviation);
    }

    @VisibleForTesting
    double maxCopyNumberDeviation(@NotNull FittedRegion fittedRegion) {
        boolean structuralBreakTransition = fittedRegion.start() < combinedRegion.region().start()
                ? !combinedRegion.region().structuralVariantSupport().equals(StructuralVariantSupport.NONE)
                : !fittedRegion.structuralVariantSupport().equals(StructuralVariantSupport.NONE);

        return maxCopyNumberDeviation(fittedRegion.bafCount(), fittedRegion.observedTumorRatioCount(), structuralBreakTransition);
    }

    private double maxCopyNumberDeviation(int bafCount, int observedTumorRatioCount, boolean structuralBreakTransition) {
        if (bafCount >= 10) {
            return MIN_COPY_NUMBER_TOLERANCE;
        }

        final double maxDeviation;
        if (structuralBreakTransition || observedTumorRatioCount > 5) {
            maxDeviation = HC_MAX_COPY_NUMBER_TOLERANCE;
        } else {
            maxDeviation = LC_MAX_COPY_NUMBER_TOLERANCE;
        }

        return (MIN_COPY_NUMBER_TOLERANCE - maxDeviation) / 10 * bafCount + maxDeviation;
    }

}
