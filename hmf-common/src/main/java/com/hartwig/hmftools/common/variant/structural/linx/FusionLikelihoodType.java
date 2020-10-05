package com.hartwig.hmftools.common.variant.structural.linx;

import org.jetbrains.annotations.NotNull;

public enum FusionLikelihoodType
{
    HIGH("High"),
    LOW("Low"),
    NA("NA");

    @NotNull
    private final String display;

    FusionLikelihoodType(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }
}
