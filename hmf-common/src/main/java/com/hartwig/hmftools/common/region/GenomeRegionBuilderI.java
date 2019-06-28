package com.hartwig.hmftools.common.region;

import org.jetbrains.annotations.NotNull;

//TODO: Rename this and GenomeRegionBuilder
public interface GenomeRegionBuilderI<T extends GenomeRegion>
{
    @NotNull
    GenomeRegionBuilderI<T> chromosome(@NotNull String chromosome);

    @NotNull
    GenomeRegionBuilderI<T> start(long start);

    @NotNull
    GenomeRegionBuilderI<T> end(long end);

    @NotNull
    T build();
}
