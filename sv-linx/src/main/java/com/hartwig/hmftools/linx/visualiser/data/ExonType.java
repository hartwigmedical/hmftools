package com.hartwig.hmftools.linx.visualiser.data;

import org.jetbrains.annotations.NotNull;

public enum ExonType
{
    DRIVER,
    FUSION,
    PSEUDOGENE,
    DISRUPTED;

    @NotNull
    public static ExonType fromString(@NotNull final String string) {
        return string.equals("EXON_LOST") ? DISRUPTED : ExonType.valueOf(string);
    }

}
