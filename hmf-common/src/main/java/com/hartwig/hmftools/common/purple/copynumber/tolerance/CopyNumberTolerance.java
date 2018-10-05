package com.hartwig.hmftools.common.purple.copynumber.tolerance;

import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

public interface CopyNumberTolerance {
    boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second);

}
