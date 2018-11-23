package com.hartwig.hmftools.amber;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.amber.NormalBAF;

public class HetrozygousFilter implements Predicate<NormalBAF> {

    private final double minHetAFPercentage;
    private final double maxHetAFPercentage;

    HetrozygousFilter(final double minHetAFPercentage, final double maxHetAFPercentage) {
        this.minHetAFPercentage = minHetAFPercentage;
        this.maxHetAFPercentage = maxHetAFPercentage;
    }

    @Override
    public boolean test(final NormalBAF bafEvidence) {
        if (bafEvidence.isValid() && isHeterozygousRef(bafEvidence.refCount(), bafEvidence.readDepth())) {
            final NormalBAF.Base alt = bafEvidence.alt();
            return isHeterozygousAlt(bafEvidence.baseMap().get(alt), bafEvidence.refCount());
        }

        return false;
    }

    private boolean isHeterozygousRef(int refCount, int totalCount) {
        final int minCount = (int) Math.round((1 - maxHetAFPercentage) * totalCount);
        final int maxCount = (int) Math.round((1 - minHetAFPercentage) * totalCount);
        return between(refCount, minCount, maxCount);
    }

    private boolean isHeterozygousAlt(int altCount, int totalCount) {
        final int minCount = (int) Math.round(minHetAFPercentage * totalCount);
        final int maxCount = (int) Math.round(maxHetAFPercentage * totalCount);
        return between(altCount, minCount, maxCount);
    }

    private static boolean between(int totalCount, int min, int max) {
        return totalCount >= min && totalCount <= max;
    }

}
