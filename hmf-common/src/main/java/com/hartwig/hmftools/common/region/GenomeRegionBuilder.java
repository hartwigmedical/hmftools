package com.hartwig.hmftools.common.region;

import org.jetbrains.annotations.NotNull;

public interface GenomeRegionBuilder<T extends GenomeRegion>
{
    @NotNull
    GenomeRegionBuilder<T> chromosome(@NotNull String chromosome);

    @NotNull
    GenomeRegionBuilder<T> start(long start);

    @NotNull
    GenomeRegionBuilder<T> end(long end);

    @NotNull
    T build();
}
