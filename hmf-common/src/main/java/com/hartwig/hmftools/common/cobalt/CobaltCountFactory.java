package com.hartwig.hmftools.common.cobalt;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

public final class CobaltCountFactory
{

    private CobaltCountFactory()
    {
    }

    @NotNull
    public static Multimap<Chromosome, CobaltCount> tumorOnly(@NotNull final Multimap<Chromosome, ReadCount> tumorCount)
    {
        final Multimap<Chromosome, CobaltCount> result = ArrayListMultimap.create();

        for(Chromosome chromosome : tumorCount.keySet())
        {
            for(ReadCount tumorReadCount : tumorCount.get(chromosome))
            {
                final CobaltCount position = fromTumor(tumorReadCount);
                result.put(chromosome, position);
            }
        }

        return result;
    }

    @NotNull
    public static Multimap<Chromosome, CobaltCount> pairedTumorNormal(@NotNull final Multimap<Chromosome, ReadCount> referenceCount,
            @NotNull final Multimap<Chromosome, ReadCount> tumorCount)
    {
        final Multimap<Chromosome, CobaltCount> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadCount> tumorCountSelector = GenomePositionSelectorFactory.create(tumorCount);

        for(Chromosome chromosome : referenceCount.keySet())
        {
            for(ReadCount referenceReadCount : referenceCount.get(chromosome))
            {
                int tumorReadCount = tumorCountSelector.select(referenceReadCount).map(ReadCount::readCount).orElse(0);
                final CobaltCount position = fromReference(referenceReadCount, tumorReadCount);
                result.put(chromosome, position);
            }
        }

        return result;
    }

    @NotNull
    private static CobaltCount fromReference(@NotNull final ReadCount reference, int tumorReadCount)
    {
        return ImmutableCobaltRatio.builder()
                .from(reference)
                .referenceReadCount(reference.readCount())
                .tumorReadCount(tumorReadCount)
                .referenceGCRatio(-1)
                .referenceGCDiploidRatio(-1)
                .tumorGCRatio(-1)
                .build();
    }

    @NotNull
    private static CobaltCount fromTumor(@NotNull final ReadCount tumor)
    {
        return ImmutableCobaltRatio.builder()
                .from(tumor)
                .tumorReadCount(tumor.readCount())
                .referenceReadCount(-1)
                .referenceGCRatio(-1)
                .referenceGCDiploidRatio(-1)
                .tumorGCRatio(-1)
                .build();
    }

}
