package com.hartwig.hmftools.common.purple.copynumber.tolerance;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

public interface CopyNumberTolerance {
    boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second);

    static CopyNumberTolerance create(boolean experimental, @NotNull final PurityAdjuster purityAdjuster) {

        if (experimental) {
            return new AlleleDeviation(purityAdjuster);
        }

        final CopyNumberDeviation copyNumberDeviation = new CopyNumberDeviation(purityAdjuster);
        return (first, second) -> copyNumberDeviation.inTolerance(first, second) && BAFDeviation.inTolerance(first, second);

    }

}
