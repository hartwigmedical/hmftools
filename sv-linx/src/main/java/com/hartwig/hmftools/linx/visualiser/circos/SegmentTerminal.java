package com.hartwig.hmftools.linx.visualiser.circos;

import org.jetbrains.annotations.NotNull;

public enum SegmentTerminal
{
    TELOMERE,
    CENTROMERE,
    NONE;

    @NotNull
    public static SegmentTerminal fromString(@NotNull final String position)
    {
        switch (position)
        {
            case "C":
                return CENTROMERE;
            case "T":
                return TELOMERE;
            default:
                return NONE;
        }
    }

}
