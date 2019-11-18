package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class ExtendStructuralVariant extends ExtendRegion {

    @NotNull
    static List<CombinedRegion> extendStructuralVariants(@NotNull final List<CombinedRegion> regions) {
        return new ExtendStructuralVariant().extend(regions);
    }

    private ExtendStructuralVariant() {
        super(CopyNumberMethod.STRUCTURAL_VARIANT);
    }

    @Override
    protected void extend(final CombinedRegion target, final CombinedRegion neighbour) {
        target.extendWithUnweightedAverage(neighbour.region());
    }
}
