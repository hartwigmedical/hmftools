package com.hartwig.hmftools.linx.visualiser.circos;

import org.jetbrains.annotations.NotNull;

public enum SegmentTerminal
{
    TELOMERE,
    CENTROMERE,
    NONE;

    public static SegmentTerminal fromString(final String position)
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
