package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;

import org.jetbrains.annotations.NotNull;

public class ExtendStructuralVariant extends ExtendRegion
{
    public static List<CombinedRegion> extendStructuralVariants(@NotNull final List<CombinedRegion> regions)
    {
        return new ExtendStructuralVariant().extend(regions);
    }

    private ExtendStructuralVariant()
    {
        super(CopyNumberMethod.STRUCTURAL_VARIANT);
    }

    @Override
    protected void extend(final CombinedRegion target, final CombinedRegion neighbour)
    {
        target.extendWithUnweightedAverage(neighbour.region());
    }
}
