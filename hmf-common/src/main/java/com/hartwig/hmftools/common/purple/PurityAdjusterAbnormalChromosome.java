package com.hartwig.hmftools.common.purple;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;

import org.jetbrains.annotations.NotNull;

public class PurityAdjusterAbnormalChromosome extends PurityAdjuster {

    private final Map<String, Double> observedRatioMap;

    public PurityAdjusterAbnormalChromosome(final double purity, final double normFactor, final List<CobaltChromosome> chromosomeList) {
        super(purity, normFactor);
        this.observedRatioMap =
                chromosomeList.stream().collect(Collectors.toMap(CobaltChromosome::contig, CobaltChromosome::observedRatio));
    }

    @Override
    public double germlineRatio(@NotNull final String contig) {
        return observedRatioMap.getOrDefault(contig, 0d);
    }
}
