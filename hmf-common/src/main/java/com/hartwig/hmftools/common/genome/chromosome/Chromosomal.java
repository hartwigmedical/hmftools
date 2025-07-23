package com.hartwig.hmftools.common.genome.chromosome;

public interface Chromosomal
{
    String chromosome();

    default HumanChromosome chr() {
        return HumanChromosome.fromString(chromosome());
    }
}
