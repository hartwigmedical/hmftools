package com.hartwig.hmftools.common.position;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public interface GenomePosition extends Comparable<GenomePosition> {

    @NotNull
    String chromosome();

    long position();

    @Override
    default int compareTo(@NotNull final GenomePosition other) {
        if (chromosome().equals(other.chromosome())) {
            return Long.compare(position(), other.position());
        }
        return HumanChromosome.fromString(chromosome()).compareTo(HumanChromosome.fromString(other.chromosome()));
    }
}
