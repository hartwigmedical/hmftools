package com.hartwig.hmftools.amber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.ImmutableAmberBAF;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

class BAFFactory {

    private final double minHetAFPercentage;
    private final double maxHetAFPercentage;
    private final double minDepthPercentage;
    private final double maxDepthPercentage;

    BAFFactory(final double minBAFPercentage, final double minHetAFPercentage, final double minDepthPercentage,
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

            if (isValidPileup(normal) && between(readCount, minDepth, maxDepth) && isHeterozygousRef(normal.referenceCount(), readCount)
                    && isHeterozygousAlt(altCount, readCount)) {
                final Character alt = alt(altCount, normal);
                tumorSelector.select(normal)
                        .filter(BAFFactory::isValidPileup)
                        .map(x -> create(alt, normal, x))
                        .filter(BAFFactory::isValidBAF)
                        .ifPresent(result::add);
            }
        }

        return result;
    }

    private static boolean isValidPileup(@NotNull final Pileup pileup) {
        return pileup.indelCount() == 0 && pileup.readCount() > 0;
    }

    private static boolean isValidBAF(@NotNull final AmberBAF baf) {
        return Doubles.isFinite(baf.tumorBAF()) & Doubles.isFinite(baf.normalBAF());
    }

    @VisibleForTesting
    static AmberBAF create(final char base, @NotNull final Pileup normal, @NotNull final Pileup tumor) {
        int tumorAltCount = tumor.mismatchCount(base);
        double tumorBaf = tumorAltCount / (double) (tumorAltCount + tumor.referenceCount());
        int normalAltCount = normal.mismatchCount(base);
        double normalBaf = normalAltCount / (double) (normalAltCount + normal.referenceCount());
        return ImmutableAmberBAF.builder()
                .from(tumor)
                .normalBAF(normalBaf)
                .normalDepth(normal.readCount())
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

    private static int maxMismatchReadCount(@NotNull final Pileup pileup) {
        return Math.max(Math.max(Math.max(pileup.gMismatchCount(), pileup.aMismatchCount()), pileup.cMismatchCount()),
                pileup.tMismatchCount());
    }

    private static boolean between(int totalCount, int min, int max) {
        return totalCount >= min && totalCount <= max;
    }

    private static int medianReadCount(@NotNull final List<Pileup> pileups) {
        final List<Integer> reads = pileups.stream().map(Pileup::readCount).sorted().collect(Collectors.toList());
        int count = reads.size();
        return count % 2 == 0 ? (reads.get(count / 2) + reads.get(count / 2 - 1)) / 2 : reads.get(count / 2);
    }
}
