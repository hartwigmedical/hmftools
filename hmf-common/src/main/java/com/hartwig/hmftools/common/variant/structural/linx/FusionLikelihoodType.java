package com.hartwig.hmftools.common.variant.structural.linx;

import org.jetbrains.annotations.NotNull;

public enum FusionLikelihoodType
{
    HIGH("High"),
    LOW("Low"),
    NA("NA");

    private final String mDisplay;

    FusionLikelihoodType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}
