package com.hartwig.hmftools.common.cobalt;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class MedianRatioFactory {

    @NotNull
    public static List<MedianRatio> create(@NotNull Multimap<Chromosome, CobaltRatio> ratios) {
        final List<MedianRatio> results = Lists.newArrayList();

        for (Chromosome contig : HumanChromosome.values()) {
            if (ratios.containsKey(contig)) {
                final List<Double> contigRatios =
                        ratios.get(contig).stream().map(CobaltRatio::referenceGCRatio).collect(Collectors.toList());
                final double medianRatio = Doubles.median(contigRatios, Doubles::positive);
                results.add(ImmutableMedianRatio.builder().chromosome(contig.toString()).medianRatio(medianRatio).build());
            }
        }
        return results;
    }

}
