package com.hartwig.hmftools.common.genome.region;

import org.jetbrains.annotations.NotNull;

public interface GenomeRegionBuilder<T extends GenomeRegion> {

    @NotNull
    GenomeRegionBuilder<T> chromosome(@NotNull String chromosome);

    @NotNull
    GenomeRegionBuilder<T> start(int start);

    @NotNull
    GenomeRegionBuilder<T> end(int end);

    @NotNull
    T build();
}
