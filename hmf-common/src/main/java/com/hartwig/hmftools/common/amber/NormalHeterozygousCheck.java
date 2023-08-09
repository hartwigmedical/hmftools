package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.utils.Doubles;

public class NormalHeterozygousCheck
{
    private final double mMinHetAFPercentage;
    private final double mMaxHetAFPercentage;

    public NormalHeterozygousCheck(final double minHetAFPercentage, final double maxHetAFPercentage)
    {
        mMinHetAFPercentage = minHetAFPercentage;
        mMaxHetAFPercentage = maxHetAFPercentage;
    }

    public boolean test(int readDepth, int refSupport, int altSupport, int indelCount)
    {
        return indelCount == 0 && altSupport > 0 && refSupport > 0
                && isHeterozygousRef(refSupport, readDepth)
                && isHeterozygousAlt(altSupport, readDepth);
    }

    private boolean isHeterozygousRef(int refSupport, int readDepth)
    {
        final double minCount = (1 - mMaxHetAFPercentage) * readDepth;
        final double maxCount = (1 - mMinHetAFPercentage) * readDepth;
        return between(refSupport, minCount, maxCount);
    }

    private boolean isHeterozygousAlt(int altSupport, int readDepth)
    {
        final double minCount = mMinHetAFPercentage * readDepth;
        final double maxCount = mMaxHetAFPercentage * readDepth;
        return between(altSupport, minCount, maxCount);
    }

    private static boolean between(int totalCount, double min, double max)
    {
        return Doubles.greaterOrEqual(totalCount, min) && Doubles.lessOrEqual(totalCount, max);
    }
}
