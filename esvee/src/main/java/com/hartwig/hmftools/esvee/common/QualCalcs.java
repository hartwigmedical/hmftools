package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.min;
import static java.lang.Math.round;

public final class QualCalcs
{
    public static int calcQual(
            final int repeatAdjustment, final int alignmentScore, final int alignmentMapQual)
    {
        double lengthFactor = (alignmentScore - repeatAdjustment)/100.0;
        return (int)round(alignmentMapQual * min(1, lengthFactor));
    }
}
