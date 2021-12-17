package com.hartwig.hmftools.common.genome.position;

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;

import org.jetbrains.annotations.NotNull;

public interface GenomePosition extends Comparable<GenomePosition> {

    @NotNull
    String chromosome();

    int position();

    @Override
    default int compareTo(@NotNull GenomePosition other) {
        int chromosomeCompare = ContigComparator.INSTANCE.compare(chromosome(), other.chromosome());
        if (chromosomeCompare == 0) {
            return Long.compare(position(), other.position());
        }

        return chromosomeCompare;
    }
}
