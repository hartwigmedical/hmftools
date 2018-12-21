package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Optional;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantLegPloidy implements StructuralVariantLeg {

    private static final double VAF_TO_USE_READ_DEPTH = 0.75;

    public abstract double observedVaf();

    public abstract double adjustedVaf();

    public abstract double weight();

    public abstract double averageImpliedPloidy();

    public abstract double unweightedImpliedPloidy();

    public abstract Optional<Double> leftCopyNumber();

    public abstract Optional<Double> rightCopyNumber();

    public double impliedRightCopyNumberWeight() {
        return canInferRight() ? weight() : 0;
    }

    public double impliedRightCopyNumber(int averageReadDepth, double averageCopyNumber) {

        if (isDecreasingFromZero(1, leftCopyNumber())) {
            return 0;
        }

        if (isIncreasingAndLargeVaf(-1)) {
            return copyNumberImpliedFromReadDepth(averageReadDepth, averageCopyNumber);
        }

        return leftCopyNumber().map(x -> x - orientation() * averageImpliedPloidy()).orElse(0D);
    }


    public double impliedLeftCopyNumber(int averageReadDepth, double averageCopyNumber) {

        if (isDecreasingFromZero(-1, rightCopyNumber())) {
            return 0;
        }

        if (isIncreasingAndLargeVaf(1)) {
            return copyNumberImpliedFromReadDepth(averageReadDepth, averageCopyNumber);
        }

        return rightCopyNumber().map(x -> x + orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedLeftCopyNumberWeight() {
        return canInferLeft() ? weight() : 0;
    }

    private boolean canInferRight() {
        return isDecreasingFromZero(1, leftCopyNumber()) || isIncreasingAndLargeVaf(-1) || leftCopyNumber().isPresent();
    }

    private boolean canInferLeft() {
        return isDecreasingFromZero(-1, rightCopyNumber()) || isIncreasingAndLargeVaf(1) || rightCopyNumber().isPresent();
    }

    private boolean isDecreasingFromZero(int orientation, @NotNull final Optional<Double> copyNumber) {
        return orientation() == orientation && copyNumber.filter(Doubles::isZero).isPresent();
    }

    private double copyNumberImpliedFromReadDepth(int averageReadDepth, double averageCopyNumber) {
        Integer tumourVariantFragmentCount = tumourVariantFragmentCount();
        Integer tumourReferenceFragmentCount = tumourReferenceFragmentCount();
        int readDepth = (tumourVariantFragmentCount == null ? 0 : tumourVariantFragmentCount) + (tumourReferenceFragmentCount == null
                ? 0
                : tumourReferenceFragmentCount);
        return averageCopyNumber * readDepth / averageReadDepth;
    }

    private boolean isIncreasingAndLargeVaf(int orientation) {
        if (Doubles.greaterOrEqual(adjustedVaf(), VAF_TO_USE_READ_DEPTH) && orientation() == orientation) {

            Integer tumourVariantFragmentCount = tumourVariantFragmentCount();
            Integer tumourReferenceFragmentCount = tumourReferenceFragmentCount();
            int readDepth = (tumourVariantFragmentCount == null ? 0 : tumourVariantFragmentCount) + (tumourReferenceFragmentCount == null
                    ? 0
                    : tumourReferenceFragmentCount);

            return readDepth > 0;
        }

        return false;
    }

}
