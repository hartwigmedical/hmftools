package com.hartwig.hmftools.common.variant.structural.linx;

import org.jetbrains.annotations.NotNull;

public enum FusionPhasedType
{
    INFRAME("Inframe"),
    SKIPPED_EXONS("Skipped exons"),
    OUT_OF_FRAME ("Out of frame");

    @NotNull
    private final String display;

    FusionPhasedType(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

}
