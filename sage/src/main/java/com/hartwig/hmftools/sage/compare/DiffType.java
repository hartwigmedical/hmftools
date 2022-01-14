package com.hartwig.hmftools.sage.compare;

import static java.lang.Math.abs;
import static java.lang.Math.max;

public enum DiffType
{
    MATCHED,
    NO_NEW,
    NO_ORIG,
    FILTER_PASS,
    FILTER_DIFF,
    QUAL,
    TIER,
    ALLELE_DEPTH,
    LOCAL_PHASE,
    OTHER_VALUE;

    public static boolean hasValueDiff(final double value1, final double value2, final double diffAbs, final double diffPerc)
    {
        if(value1 == 0 && value2 == 0)
            return false;

        double diff = abs(value1 - value2);
        return diff > diffAbs && diff / max(value1, value2) > diffPerc;
    }

}
