package com.hartwig.hmftools.common.amber;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.utils.Doubles;

public class NormalHeterozygousFilter implements Predicate<BaseDepth>
{
    private final double mMinHetAFPercentage;
    private final double mMaxHetAFPercentage;

    public NormalHeterozygousFilter(final double minHetAFPercentage, final double maxHetAFPercentage)
    {
        mMinHetAFPercentage = minHetAFPercentage;
        mMaxHetAFPercentage = maxHetAFPercentage;
    }

    @Override
    public boolean test(final BaseDepth bafEvidence)
    {
        return bafEvidence.isValid() && bafEvidence.AltSupport > 0 && bafEvidence.RefSupport > 0
                && isHeterozygousRef(bafEvidence.RefSupport, bafEvidence.ReadDepth)
                && isHeterozygousAlt(bafEvidence.AltSupport, bafEvidence.ReadDepth);
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
