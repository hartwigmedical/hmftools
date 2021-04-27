package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Optional;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantLegPloidy implements StructuralVariantLeg, StructuralVariantLegCopyNumber
{
    public abstract double observedVaf();

    public abstract double adjustedVaf();

    public abstract double weight();

    public abstract double averageImpliedPloidy();

    public abstract double unweightedImpliedPloidy();

    public double impliedRightCopyNumberWeight()
    {
        return canInferRight() ? weight() : 0;
    }

    public double impliedRightCopyNumber()
    {
        if(isDecreasingFromZero(1, leftCopyNumber()))
            return 0;

        return leftCopyNumber().map(x -> x - orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedLeftCopyNumber()
    {
        if(isDecreasingFromZero(-1, rightCopyNumber()))
            return 0;

        return rightCopyNumber().map(x -> x + orientation() * averageImpliedPloidy()).orElse(0D);
    }

    public double impliedLeftCopyNumberWeight()
    {
        return canInferLeft() ? weight() : 0;
    }

    private boolean canInferRight()
    {
        return isDecreasingFromZero(1, leftCopyNumber()) || leftCopyNumber().isPresent();
    }

    private boolean canInferLeft()
    {
        return isDecreasingFromZero(-1, rightCopyNumber()) || rightCopyNumber().isPresent();
    }

    private boolean isDecreasingFromZero(int orientation, @NotNull final Optional<Double> copyNumber)
    {
        return orientation() == orientation && copyNumber.filter(Doubles::isZero).isPresent();
    }
}
