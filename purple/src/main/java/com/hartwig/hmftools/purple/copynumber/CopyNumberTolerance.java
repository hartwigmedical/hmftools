package com.hartwig.hmftools.purple.copynumber;

import com.hartwig.hmftools.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

public interface CopyNumberTolerance {

    boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second);

}
