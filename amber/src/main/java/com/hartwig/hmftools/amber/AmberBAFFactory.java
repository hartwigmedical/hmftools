package com.hartwig.hmftools.amber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.ImmutableAmberBAF;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

class AmberBAFFactory {

    private final double minHetAFPercentage;
    private final double maxHetAFPercentage;
    private final double minDepthPercentage;
    private final double maxDepthPercentage;

    AmberBAFFactory(final double minBAFPercentage, final double minHetAFPercentage, final double minDepthPercentage,
            final double maxDepthPercentage) {
        this.minHetAFPercentage = minBAFPercentage;
        this.maxHetAFPercentage = minHetAFPercentage;
        this.minDepthPercentage = minDepthPercentage;
        this.maxDepthPercentage = maxDepthPercentage;
    }

    @NotNull
    List<AmberBAF> create(@NotNull final List<Pileup> normalPileup, @NotNull final List<Pileup> tumorPileup) {

        final GenomePositionSelector<Pileup> tumorSelector = GenomePositionSelectorFactory.create(tumorPileup);

        int medianDepth = medianReadCount(normalPileup);
        int minDepth = (int) Math.round(medianDepth * minDepthPercentage);
        int maxDepth = (int) Math.round(medianDepth * maxDepthPercentage);

        final List<AmberBAF> result = Lists.newArrayList();
        for (Pileup normal : normalPileup) {

            final int readCount = normal.readCount();
            final int altCount = maxMismatchReadCount(normal);

            if (normal.indels() == 0 && between(readCount, minDepth, maxDepth) && isHetrozygousRef(normal.referenceCount(), readCount)
                    && isHetrozygousAlt(altCount, readCount)) {

                final Character alt = alt(altCount, normal);
                tumorSelector.select(normal).filter(x -> x.indels() == 0).map(x -> create(alt, normal, x)).ifPresent(result::add);
            }
        }

        return result;

    }

    private static AmberBAF create(final char base, @NotNull final Pileup normal, @NotNull final Pileup tumor) {
        int tumorAltCount = tumor.mismatchCount(base);
        double tumorBaf = tumorAltCount / (double) (tumorAltCount + tumor.referenceCount());
        int normalAltCount = normal.mismatchCount(base);
        double normalBaf = normalAltCount / (double) (normalAltCount + normal.referenceCount());
        return ImmutableAmberBAF.builder()
                .from(tumor)
                .normalBAF(normalBaf)
                .normalDepth(normal.referenceCount())
                .tumorBAF(tumorBaf)
                .tumorDepth(tumor.readCount())
                .build();
    }

    private static char alt(int count, @NotNull final Pileup pileup) {
        if (count == pileup.gMismatchCount()) {
            return 'G';
        }
        if (count == pileup.aMismatchCount()) {
            return 'A';
        }
        if (count == pileup.tMismatchCount()) {
            return 'T';
        }
        if (count == pileup.cMismatchCount()) {
            return 'C';
        }

        throw new IllegalArgumentException("unable to find count");
    }

    private boolean isHetrozygousRef(int refCount, int totalCount) {
        final int minCount = (int) Math.round((1 - maxHetAFPercentage) * totalCount);
        final int maxCount = (int) Math.round((1 - minHetAFPercentage) * totalCount);
        return between(refCount, minCount, maxCount);
    }

    private boolean isHetrozygousAlt(int altCount, int totalCount) {
        final int minCount = (int) Math.round(minHetAFPercentage * totalCount);
        final int maxCount = (int) Math.round(maxHetAFPercentage * totalCount);
        return between(altCount, minCount, maxCount);
    }

    private static int maxMismatchReadCount(@NotNull final Pileup pileup) {
        return Math.max(Math.max(Math.max(pileup.gMismatchCount(), pileup.aMismatchCount()), pileup.cMismatchCount()),
                pileup.tMismatchCount());
    }

    private boolean between(int totalCount, int min, int max) {
        return totalCount >= min && totalCount <= max;
    }

    private int medianReadCount(@NotNull final List<Pileup> pileups) {
        final List<Integer> reads = pileups.stream().map(Pileup::readCount).sorted().collect(Collectors.toList());
        int count = reads.size();
        return count % 2 == 0 ? (reads.get(count / 2) + reads.get(count / 2 - 1)) / 2 : reads.get(count / 2);
    }
}
