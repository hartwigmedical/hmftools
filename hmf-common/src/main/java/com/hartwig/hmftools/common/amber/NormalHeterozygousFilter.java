package com.hartwig.hmftools.common.amber;

import java.util.function.Predicate;

public class NormalHeterozygousFilter implements Predicate<BaseDepth>
{
    private final double minHetAFPercentage;
    private final double maxHetAFPercentage;

    public NormalHeterozygousFilter(final double minHetAFPercentage, final double maxHetAFPercentage)
    {
        this.minHetAFPercentage = minHetAFPercentage;
        this.maxHetAFPercentage = maxHetAFPercentage;
    }

    @Override
    public boolean test(final BaseDepth bafEvidence)
    {
        return bafEvidence.isValid() && bafEvidence.altSupport() > 0 && bafEvidence.refSupport() > 0
                && isHeterozygousRef(bafEvidence.refSupport(), bafEvidence.readDepth())
                && isHeterozygousAlt(bafEvidence.altSupport(), bafEvidence.readDepth());
    }

    private boolean isHeterozygousRef(int refSupport, int readDepth)
    {
        final int minCount = (int) Math.round((1 - maxHetAFPercentage) * readDepth);
        final int maxCount = (int) Math.round((1 - minHetAFPercentage) * readDepth);
        return between(refSupport, minCount, maxCount);
    }

    private boolean isHeterozygousAlt(int altSupport, int readDepth)
    {
        final int minCount = (int) Math.round(minHetAFPercentage * readDepth);
        final int maxCount = (int) Math.round(maxHetAFPercentage * readDepth);
        return between(altSupport, minCount, maxCount);
    }

    private static boolean between(int totalCount, int min, int max)
    {
        return totalCount >= min && totalCount <= max;
    }
}
