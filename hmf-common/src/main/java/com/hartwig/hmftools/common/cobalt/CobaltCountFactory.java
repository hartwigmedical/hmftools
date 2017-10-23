package com.hartwig.hmftools.common.cobalt;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CobaltCountFactory {

    private static final Logger LOGGER = LogManager.getLogger(CobaltCountFactory.class);

    @NotNull
    public static Multimap<Chromosome, CobaltCount> merge(@NotNull final Multimap<String, ReadCount> referenceCount,
            @NotNull final Multimap<String, ReadCount> tumorCount) {
        final Multimap<Chromosome, CobaltCount> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadCount> tumorCountSelector = GenomePositionSelectorFactory.create(tumorCount);

        for (String chromosomeName : referenceCount.keySet()) {
            if (HumanChromosome.contains(chromosomeName)) {
                Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
                for (ReadCount referenceReadCount : referenceCount.get(chromosomeName)) {
                    int tumorReadCount = tumorCountSelector.select(referenceReadCount).map(ReadCount::readCount).orElse(0);
                    final CobaltCount position = create(referenceReadCount, tumorReadCount);
                    result.put(chromosome, position);
                }
            } else {
                LOGGER.info("Excluding unsupported {} chromosome from read counts", chromosomeName);

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
                .tumorGCRatio(-1).build();
    }
}
