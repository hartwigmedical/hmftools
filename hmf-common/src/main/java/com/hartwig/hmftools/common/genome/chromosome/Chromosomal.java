package com.hartwig.hmftools.common.genome.chromosome;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public interface Chromosomal
{
    String chromosome();

    default HumanChromosome chr() {
        return HumanChromosome.fromString(chromosome());
    }

    static <T extends Chromosomal> Map<Chromosome, List<T>> toPerChromosomeLists(List<T> data)
    {
        final Map<Chromosome, List<T>> chrRegionsMap = Maps.newHashMap();

        for (T t : data) {
            chrRegionsMap.computeIfAbsent(t.chr(), k -> new ArrayList<>());
            chrRegionsMap.get(t.chr()).add(t);
        }
        return chrRegionsMap;
    }
}
