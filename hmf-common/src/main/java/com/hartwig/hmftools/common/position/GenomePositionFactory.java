package com.hartwig.hmftools.common.position;

import org.jetbrains.annotations.NotNull;

public class GenomePositionFactory {

    private GenomePositionFactory() {
    }

    @NotNull
    public static GenomePosition create(@NotNull final String chromosome, final long position) {
        return ImmutableGenomePositionImpl.builder().chromosome(chromosome).position(position).build();
    }

}
