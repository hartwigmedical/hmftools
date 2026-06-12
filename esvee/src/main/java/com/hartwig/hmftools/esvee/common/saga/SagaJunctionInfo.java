package com.hartwig.hmftools.esvee.common.saga;

import com.hartwig.hmftools.common.genome.region.Orientation;

import org.jetbrains.annotations.Nullable;

public record SagaJunctionInfo(
        int assemblyOffset,
        @Nullable Orientation orientation
)
{
}
