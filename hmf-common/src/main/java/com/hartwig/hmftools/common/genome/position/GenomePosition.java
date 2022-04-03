package com.hartwig.hmftools.common.genome.position;

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;

import org.jetbrains.annotations.NotNull;

public interface GenomePosition extends Comparable<GenomePosition> {

    @NotNull
    String chromosome();

    int position();

    @Override
    default int compareTo(@NotNull GenomePosition other) {
        return compare(this, other);
    }

    // this makes it easier to use the following compare function in Map, for classes that
    // extends GenomePosition but may have their own compareTo function.
    // i.e.
    // new TreeMap(GenomePosition::compare);
    static int compare(final GenomePosition gp1, final GenomePosition gp2) {
        int chromosomeCompare = ContigComparator.INSTANCE.compare(gp1.chromosome(), gp2.chromosome());
        if (chromosomeCompare == 0) {
            return Integer.compare(gp1.position(), gp2.position());
        }
        return chromosomeCompare;
    }
}
