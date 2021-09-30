package com.hartwig.hmftools.common.purple;

import java.util.Collection;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;

import org.jetbrains.annotations.NotNull;

public class PurityAdjusterAbnormalChromosome extends PurityAdjuster
{
    private final Map<String, Double> mObservedRatioMap;

    public PurityAdjusterAbnormalChromosome(final double purity, final double normFactor,
            final Collection<CobaltChromosome> chromosomeList)
    {
        super(purity, normFactor);
        mObservedRatioMap = chromosomeList.stream().collect(Collectors.toMap(CobaltChromosome::contig, CobaltChromosome::actualRatio));
    }

    @Override
    public double germlineRatio(@NotNull final String contig)
    {
        return mObservedRatioMap.getOrDefault(contig, 0d);
    }
}
