package com.hartwig.hmftools.common.cobalt;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

public final class CobaltCountFactory {

    @NotNull
    public static Multimap<Chromosome, CobaltCount> merge(@NotNull final Multimap<Chromosome, ReadCount> referenceCount,
            @NotNull final Multimap<Chromosome, ReadCount> tumorCount) {
        final Multimap<Chromosome, CobaltCount> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadCount> tumorCountSelector = GenomePositionSelectorFactory.create(tumorCount);

        for (Chromosome chromosome : referenceCount.keySet()) {
            for (ReadCount referenceReadCount : referenceCount.get(chromosome)) {
                int tumorReadCount = tumorCountSelector.select(referenceReadCount).map(ReadCount::readCount).orElse(0);
                final CobaltCount position = create(referenceReadCount, tumorReadCount);
                result.put(chromosome, position);
            }
        }

        return result;
    }

    @NotNull
    private static CobaltCount create(@NotNull final ReadCount reference, int tumorReadCount) {
        return ImmutableCobaltRatio.builder()
                .from(reference)
                .referenceReadCount(reference.readCount())
                .tumorReadCount(tumorReadCount)
                .referenceGCRatio(-1)
                .referenceGCDiploidRatio(-1)
                .tumorGCRatio(-1)
                .build();
    }
}
