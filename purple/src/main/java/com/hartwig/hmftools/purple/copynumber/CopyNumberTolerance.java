package com.hartwig.hmftools.purple.copynumber;

import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public interface CopyNumberTolerance {

    boolean inTolerance(@NotNull final ObservedRegion first, @NotNull final ObservedRegion second);

}
