package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;

abstract class BaseCopyNumberBuilder {

    private static final double MIN_COPY_NUMBER_TOLERANCE = 0.3;
    @VisibleForTesting
    static final double LC_MAX_COPY_NUMBER_TOLERANCE = 1.3;
    @VisibleForTesting
    static final double HC_MAX_COPY_NUMBER_TOLERANCE = 0.7;

    @NotNull
    private final PurityAdjuster purityAdjuster;
    private final String chromosome;
    private long start = 1;
    private long end;
    private boolean ratioSupport;
    private StructuralVariantSupport structuralVariantSupport;

    BaseCopyNumberBuilder(@NotNull final PurityAdjuster purityAdjuster, @NotNull final FittedRegion fittedRegion) {
        this.purityAdjuster = purityAdjuster;
        this.chromosome = fittedRegion.chromosome();
        this.start = fittedRegion.start();
        this.structuralVariantSupport = fittedRegion.structuralVariantSupport();
        this.ratioSupport = fittedRegion.ratioSupport();
        extendRegion(fittedRegion);
    }

    public String chromosome() {
        return chromosome;
    }

    public abstract int bafCount();

    public abstract double averageObservedBAF();

    public abstract double averageTumorCopyNumber();

    abstract double averageRefNormalisedCopyNumber();

    public void extendRegion(@NotNull final FittedRegion value) {
        assert (chromosome.equals(value.chromosome())) : "Regions cannot be extended between chromosomes";

        if (value.start() <= start) {
            structuralVariantSupport = value.structuralVariantSupport();
            ratioSupport = value.ratioSupport();
        }

        start = Math.min(value.start(), start);
        end = Math.max(value.end(), end);
    }

    @NotNull
    public PurpleCopyNumber build() {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .bafCount(bafCount())
                .averageObservedBAF(averageObservedBAF())
                .averageActualBAF(purityAdjustedBAF(averageObservedBAF()))
                .averageTumorCopyNumber(averageTumorCopyNumber())
                .ratioSupport(ratioSupport)
                .structuralVariantSupport(structuralVariantSupport)
                .build();
    }

    @VisibleForTesting
    double purityAdjustedBAF(final double observedBAF) {

        double copyNumber = averageTumorCopyNumber();

        if (Doubles.isZero(copyNumber)) {
            return observedBAF;
        }

        return purityAdjuster.purityAdjustedBAF(chromosome, copyNumber, observedBAF);
    }

    boolean withinCopyNumberTolerance(@NotNull final FittedRegion copyNumber) {
        double tumorCopyNumberDeviation = Math.abs(copyNumber.tumorCopyNumber() - averageTumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(copyNumber.refNormalisedCopyNumber() - averageRefNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);
        return Doubles.lessOrEqual(copyNumberDeviation, maxCopyNumberDeviation(copyNumber));
    }

    @VisibleForTesting
    double maxCopyNumberDeviation(@NotNull FittedRegion fittedRegion) {
        if (fittedRegion.bafCount() >= 10) {
            return MIN_COPY_NUMBER_TOLERANCE;
        }

        boolean structuralBreakTransition = fittedRegion.start() < start
                ? !structuralVariantSupport.equals(StructuralVariantSupport.NONE)
                : !fittedRegion.structuralVariantSupport().equals(StructuralVariantSupport.NONE);

        final double maxDeviation;
        if (structuralBreakTransition || fittedRegion.observedTumorRatioCount() > 5) {
            maxDeviation = HC_MAX_COPY_NUMBER_TOLERANCE;
        } else {
            maxDeviation = LC_MAX_COPY_NUMBER_TOLERANCE;
        }

        double result = (MIN_COPY_NUMBER_TOLERANCE - maxDeviation) / 10 * fittedRegion.bafCount() + maxDeviation;

        // Adjust for low purity
        return result * Math.max(1, 0.20 / purityAdjuster.purity());
    }

}
