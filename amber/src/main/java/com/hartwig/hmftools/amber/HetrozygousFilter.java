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
        return bafEvidence.isValid()
                && isHeterozygousRef(bafEvidence.refSupport(), bafEvidence.readDepth())
                && isHeterozygousAlt(bafEvidence.altSupport(), bafEvidence.readDepth());
    }

    private boolean isHeterozygousRef(int refSupport, int readDepth) {
        final int minCount = (int) Math.round((1 - maxHetAFPercentage) * readDepth);
        final int maxCount = (int) Math.round((1 - minHetAFPercentage) * readDepth);
        return between(refSupport, minCount, maxCount);
    }

    private boolean isHeterozygousAlt(int altSupport, int readDepth) {
        final int minCount = (int) Math.round(minHetAFPercentage * readDepth);
        final int maxCount = (int) Math.round(maxHetAFPercentage * readDepth);
        return between(altSupport, minCount, maxCount);
    }

    private static boolean between(int totalCount, int min, int max) {
        return totalCount >= min && totalCount <= max;
    }

}
