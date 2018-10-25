package com.hartwig.hmftools.common.collect;

import java.util.Collection;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class Multimaps {

    @NotNull
    public static <T extends GenomeRegion> ListMultimap<Chromosome, T> fromRegions(@NotNull final Collection<T> regions) {
        final ListMultimap<Chromosome, T> result = ArrayListMultimap.create();
        for (T region : regions) {
            if (HumanChromosome.contains(region.chromosome())) {
                result.put(HumanChromosome.fromString(region.chromosome()), region);
            }
        }

        return result;
    }

    @NotNull
    public static <T extends GenomePosition> ListMultimap<Chromosome, T> fromPositions(@NotNull final Collection<T> regions) {
        final ListMultimap<Chromosome, T> result = ArrayListMultimap.create();
        for (T region : regions) {
            if (HumanChromosome.contains(region.chromosome())) {
                result.put(HumanChromosome.fromString(region.chromosome()), region);
            }
        }

        return result;
    }
}
