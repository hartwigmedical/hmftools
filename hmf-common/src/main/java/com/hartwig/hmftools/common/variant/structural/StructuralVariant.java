package com.hartwig.hmftools.common.variant.structural;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface StructuralVariant {

    @Nullable
    Integer primaryKey();

    @NotNull
    String id();

    @Nullable
    String mateId();

    @NotNull
    String startChromosome();

    @NotNull
    String endChromosome();

    default String chromosome(final boolean isStart) {
        return isStart ? startChromosome() : endChromosome();
    }

    long startPosition();

    long endPosition();

    default long position(final boolean isStart) {
        return isStart ? startPosition() : endPosition();
    }

    byte startOrientation();

    byte endOrientation();

    @NotNull
    String startHomology();

    @NotNull
    String endHomology();

    @Nullable
    Double startAF();

    @Nullable
    Double endAF();

    @NotNull
    String insertSequence();

    @NotNull
    StructuralVariantType type();
}
