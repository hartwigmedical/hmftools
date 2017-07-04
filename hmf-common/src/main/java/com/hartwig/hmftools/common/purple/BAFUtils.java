package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.numeric.Doubles;

public class BAFUtils {

    public static final double NORMAL_BAF = 0.535;

    public static int minAlleleCount(final int ploidy) {
        return (int) Math.max(0, Math.round(ploidy / 2d));
    }

    public static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        assert (alleleCount >= ploidy / 2d);

        double divisor =  (2 + purity * (ploidy - 2));
        return Doubles.isZero(divisor) ? NORMAL_BAF : Math.max(NORMAL_BAF, (1 + purity * (alleleCount - 1)) / divisor);
    }

}
