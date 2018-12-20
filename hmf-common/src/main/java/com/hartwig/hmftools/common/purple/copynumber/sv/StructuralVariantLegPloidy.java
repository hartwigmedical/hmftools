package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Optional;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantLegPloidy implements StructuralVariantLeg {

    private static final double VAF_TO_USE_READ_DEPTH = 0.75;

    public abstract double vaf();

    public abstract Optional<Double> leftCopyNumber();

    public double impliedRightCopyNumber(int averageReadDepth, double averageCopyNumber) {

        if (useReadDepth(-1)) {
            return copyNumberImpliedFromReadDepth(averageReadDepth, averageCopyNumber);
        }

        return leftCopyNumber().map(x -> x - orientation() * averageImpliedPloidy()).orElse(0D);
    }

    @Deprecated
    public double impliedRightCopyNumber() {
        return leftCopyNumber().map(x -> x - orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedRightCopyNumberWeight() {
        return useReadDepth(-1) || leftCopyNumber().isPresent() ? weight() : 0;
    }

    public abstract Optional<Double> rightCopyNumber();

    public double impliedLeftCopyNumber(int averageReadDepth, double averageCopyNumber) {
        if (useReadDepth(1)) {
            return copyNumberImpliedFromReadDepth(averageReadDepth, averageCopyNumber);
        }

        return rightCopyNumber().map(x -> x + orientation() * averageImpliedPloidy()).orElse(0D);
    }

    @Deprecated
    public double impliedLeftCopyNumber() {
        return rightCopyNumber().map(x -> x + orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedLeftCopyNumberWeight() {
        return useReadDepth(1) || rightCopyNumber().isPresent() ? weight() : 0;
    }

    public abstract double unweightedImpliedPloidy();

    public abstract double averageImpliedPloidy();

    public abstract double weight();

    private double copyNumberImpliedFromReadDepth(int averageReadDepth, double averageCopyNumber) {
        Integer tumourVariantFragmentCount = tumourVariantFragmentCount();
        Integer tumourReferenceFragmentCount = tumourReferenceFragmentCount();
        int readDepth = (tumourVariantFragmentCount == null ? 0 : tumourVariantFragmentCount) + (tumourReferenceFragmentCount == null
                ? 0
                : tumourReferenceFragmentCount);
        return averageCopyNumber * readDepth / averageReadDepth;
    }

    private boolean useReadDepth(int orientation) {
        if (Doubles.greaterOrEqual(vaf(), VAF_TO_USE_READ_DEPTH) && orientation() == orientation) {

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
