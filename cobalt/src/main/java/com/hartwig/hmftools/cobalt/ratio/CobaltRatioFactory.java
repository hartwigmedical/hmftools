package com.hartwig.hmftools.cobalt.ratio;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

public final class CobaltRatioFactory
{

    private CobaltRatioFactory()
    {
    }

    @NotNull
    public static Multimap<Chromosome, CobaltRatio> merge(@NotNull final Multimap<Chromosome, CobaltCount> counts,
            @NotNull final Multimap<Chromosome, ReadRatio> referenceGCRatio, @NotNull final Multimap<Chromosome, ReadRatio> tumorGCRatio,
            @NotNull final Multimap<Chromosome, ReadRatio> referenceGCDiploidRatio)
    {
        final Multimap<Chromosome, CobaltRatio> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadRatio> tumorGCRatioSelector = GenomePositionSelectorFactory.create(tumorGCRatio);
        final GenomePositionSelector<ReadRatio> referenceGCRatioSelector = GenomePositionSelectorFactory.create(referenceGCRatio);
        final GenomePositionSelector<ReadRatio> referenceGCDiploidRatioSelector =
                GenomePositionSelectorFactory.create(referenceGCDiploidRatio);
        for(Chromosome chromosome : counts.keySet())
        {
            for(CobaltCount count : counts.get(chromosome))
            {
                final CobaltRatio ratio = ImmutableCobaltRatio.builder()
                        .from(count)
                        .tumorGCRatio(tumorGCRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .referenceGCRatio(referenceGCRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .referenceGCDiploidRatio(referenceGCDiploidRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .build();

                result.put(chromosome, ratio);
            }
        }

        return result;
    }
}
