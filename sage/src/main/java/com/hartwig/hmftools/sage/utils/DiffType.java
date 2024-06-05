package com.hartwig.hmftools.sage.utils;

import static java.lang.Math.abs;
import static java.lang.Math.max;

public enum DiffType
{
    NO_NEW,
    NO_ORIG,
    VALUE;

    public static boolean hasValueDiff(final double value1, final double value2, final double diffAbs, final double diffPerc)
    {
        if(value1 == 0 && value2 == 0)
            return false;

        double diff = abs(value1 - value2);
        return diff > diffAbs && diff / max(value1, value2) > diffPerc;
    }

}
