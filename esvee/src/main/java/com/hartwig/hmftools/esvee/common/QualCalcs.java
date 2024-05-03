package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.min;
import static java.lang.Math.round;

public final class QualCalcs
{
    public static int calcQual()
    {
        /*
        int repeatAdjustment = Alignment.segmentLength() - Alignment.repeatTrimmedLength();
        double lengthFactor = (Alignment.Score - repeatAdjustment)/100.0;
        return (int)round(Alignment.MapQual * min(1, lengthFactor));
         */

        return 0;
    }

}
