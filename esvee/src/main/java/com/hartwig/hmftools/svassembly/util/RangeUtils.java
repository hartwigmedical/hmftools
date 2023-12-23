package com.hartwig.hmftools.svassembly.util;

public enum RangeUtils
{
    ;

    public static boolean touches(final int leftStart, final int leftEnd, final int rightStart, final int rightEnd)
    {
        return leftStart - 1 <= rightEnd && leftEnd + 1 >= rightStart;
    }

    public static boolean overlaps(final int leftStart, final int leftEnd, final int rightStart, final int rightEnd)
    {
        return leftStart <= rightEnd && leftEnd >= rightStart;
    }
}
