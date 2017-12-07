package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantLegPloidy implements StructuralVariantLeg {

    public abstract double vaf();

    public abstract Optional<Double> leftCopyNumber();

    public double impliedRightCopyNumber() {
        return leftCopyNumber().map(x -> x - orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedRightCopyNumberWeight() {
        return leftCopyNumber().isPresent() ? weight() : 0;
    }

    public abstract Optional<Double> rightCopyNumber();

    public double impliedLeftCopyNumber() {
        return rightCopyNumber().map(x -> x + orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedLeftCopyNumberWeight() {
        return rightCopyNumber().isPresent() ? weight() : 0;
    }

    public abstract double unweightedImpliedPloidy();

    public abstract double averageImpliedPloidy();

    public abstract double weight();
}
