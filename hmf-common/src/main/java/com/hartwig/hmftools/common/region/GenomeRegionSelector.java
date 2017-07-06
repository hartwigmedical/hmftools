package com.hartwig.hmftools.common.region;

import java.util.Optional;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public interface GenomeRegionSelector<R extends GenomeRegion> {
    @NotNull
    Optional<R> select(GenomePosition position);
}
