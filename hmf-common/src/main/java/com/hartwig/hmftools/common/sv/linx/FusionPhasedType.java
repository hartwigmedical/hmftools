package com.hartwig.hmftools.common.sv.linx;

import org.jetbrains.annotations.NotNull;

public enum FusionPhasedType
{
    INFRAME("Inframe"),
    SKIPPED_EXONS("Skipped exons"),
    OUT_OF_FRAME ("Out of frame");

    private final String mDisplay;

    FusionPhasedType(@NotNull final String display)
    {
        mDisplay = display;
    }

    public String display()
    {
        return mDisplay;
    }

}
