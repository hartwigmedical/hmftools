package com.hartwig.hmftools.common.position;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class GenomePositionSelectorFactory {
    @NotNull
    public static <P extends GenomePosition> GenomePositionSelector<P> create(@NotNull final List<P> positions) {
        return new GenomePositionSelectorImpl<>(positions);
    }

}
