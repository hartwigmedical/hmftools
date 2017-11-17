package com.hartwig.hmftools.common.region;

import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public interface GenomeRegionSelector<R extends GenomeRegion> {
    @NotNull
    Optional<R> select(@NotNull GenomePosition position);

    void select(@NotNull GenomeRegion region, @NotNull Consumer<R> handler);
}
