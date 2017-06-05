package com.hartwig.hmftools.common.position;

import com.hartwig.hmftools.common.chromosome.Chromosomes;

import org.jetbrains.annotations.NotNull;

public interface GenomePosition extends Comparable<GenomePosition> {

    @NotNull
    String chromosome();

    long position();

    @NotNull
    default String chromosomePosition() {
        return chromosome() + ":" + position();
    }

    @Override
    default int compareTo(@NotNull final GenomePosition other) {
        if (chromosome().equals(other.chromosome())) {
            return Long.compare(position(), other.position());
        }

        return Integer.compare(Chromosomes.asInt(chromosome()), Chromosomes.asInt(other.chromosome()));
    }
}
