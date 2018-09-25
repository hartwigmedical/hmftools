package com.hartwig.hmftools.common.purple.copynumber.tolerance;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

public class BAFDeviation {

    public static boolean inTolerance(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        int minBafCount = Math.min(left.bafCount(), right.bafCount());

        if (left.bafCount() > 0 && right.bafCount() > 0) {
            double bafDeviation = Math.abs(left.observedBAF() - right.observedBAF());
            if (Doubles.greaterThan(bafDeviation, allowedBAFDeviation(minBafCount))) {
                return false;
            }
        }

        return true;
    }

    private static double allowedBAFDeviation(int bafCount) {
        return 1d / Math.max(1, bafCount) * 0.35 + 0.03;
    }
}
