package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purity.PurityAdjustment;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;

abstract class BaseCopyNumberBuilder {

    private static final double MIN_COPY_NUMBER_TOLERANCE = 0.3;
    @VisibleForTesting
    static final double MAX_COPY_NUMBER_TOLERANCE = 1.3;
    @VisibleForTesting
    static final double STRUCTURAL_VARIANCE_MAX_COPY_NUMBER_TOLERANCE = 0.7;

    private final double purity;
    private final String chromosome;
    private long start = 1;
    private long end;
    private boolean ratioSupport;
    private StructuralVariantSupport structuralVariantSupport;

    BaseCopyNumberBuilder(double purity, @NotNull final FittedRegion fittedRegion) {
        this.purity = purity;
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
                .averageActualBAF(purityAdjustedBAF(purity, averageTumorCopyNumber(), averageObservedBAF()))
                .averageTumorCopyNumber(averageTumorCopyNumber())
                .ratioSupport(ratioSupport)
                .structuralVariantSupport(structuralVariantSupport)
                .build();
    }

    @VisibleForTesting
    static double purityAdjustedBAF(final double purity, final double copyNumber, final double observedBAF) {
        if (Doubles.isZero(copyNumber)) {
            return observedBAF;
        }

        double adjustedObservedBAF = isEven(copyNumber) && lessOrEqual(observedBAF, FittedRegionFactory.NORMAL_BAF) ? 0.5 : observedBAF;
        return PurityAdjustment.purityAdjustedBAF(purity, copyNumber, adjustedObservedBAF);
    }

    @VisibleForTesting
    static boolean isEven(double copyNumber) {

        double decimal = copyNumber % 1d;
        double wholeNumber = copyNumber - decimal;

        return (wholeNumber % 2 == 0 && Doubles.lessOrEqual(decimal, 0.25)) || (wholeNumber % 2 != 0 && Doubles.greaterOrEqual(decimal,
                0.75));
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

        double maxDeviation = structuralBreakTransition ? STRUCTURAL_VARIANCE_MAX_COPY_NUMBER_TOLERANCE : MAX_COPY_NUMBER_TOLERANCE;
        return (MIN_COPY_NUMBER_TOLERANCE - maxDeviation) / 10 * fittedRegion.bafCount() + maxDeviation;
    }

}
