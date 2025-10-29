package com.hartwig.hmftools.common.genome.position;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.Chromosomal;
import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface GenomePosition extends Chromosomal, Comparable<GenomePosition>
{
    int position();

    @Override
    default int compareTo(@NotNull GenomePosition other)
    {
        return compare(this, other);
    }

    default  <T extends GenomeRegion> List<T> findContainingRegions(@NotNull List<T> regions)
    {
        return regions.stream().filter(region -> region.contains(this)).collect(Collectors.toList());
    }

    // this makes it easier to use the following compare function in Map, for classes that
    // extends GenomePosition but may have their own compareTo function.
    // i.e.
    // new TreeMap(GenomePosition::compare);
    static int compare(final GenomePosition gp1, final GenomePosition gp2)
    {
        int chromosomeCompare = ContigComparator.INSTANCE.compare(gp1.chromosome(), gp2.chromosome());
        if(chromosomeCompare == 0)
        {
            return Integer.compare(gp1.position(), gp2.position());
        }
        return chromosomeCompare;
    }
}
