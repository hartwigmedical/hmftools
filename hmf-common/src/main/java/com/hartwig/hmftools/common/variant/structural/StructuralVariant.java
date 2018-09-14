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
    StructuralVariantLeg start();

    @Nullable
    StructuralVariantLeg end();

    default String chromosome(final boolean isStart) {
        return isStart ? start().chromosome() : end().chromosome();
    }

    default long position(final boolean isStart) {
        return isStart ? start().position() : end().position();
    }

    default byte orientation(final boolean isStart) {
        return  isStart ? start().orientation() : end().orientation();
    }

    @NotNull
    String insertSequence();

    @NotNull
    StructuralVariantType type();

    @Nullable
    String filter();

    @Nullable
    Boolean imprecise();

    @Nullable
    Double qualityScore();

    @Nullable
    String event();

    @Nullable
    String startLinkedBy();

    @Nullable
    String endLinkedBy();
}
