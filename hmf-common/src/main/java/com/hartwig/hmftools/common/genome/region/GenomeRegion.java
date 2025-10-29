package com.hartwig.hmftools.common.genome.region;

import com.hartwig.hmftools.common.genome.chromosome.Chromosomal;
import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegionProvider;

public interface GenomeRegion extends Chromosomal, Comparable<GenomeRegion>, ChrBaseRegionProvider
{
    int start();

    int end();

    default int bases()
    {
        return 1 + end() - start();
    }

    @Override
    default int compareTo(final GenomeRegion other)
    {
        if(chromosome().equals(other.chromosome()))
        {
            if(start() < other.start())
            {
                return -1;
            }
            else if(start() == other.start())
            {
                return 0;
            }
            return 1;
        }

        return ContigComparator.INSTANCE.compare(chromosome(), other.chromosome());
    }

    default boolean contains(final GenomePosition position)
    {
        return chromosome().equals(position.chromosome()) && start() <= position.position() && end() >= position.position();
    }

    default boolean overlaps(final GenomeRegion other)
    {
        return other.chromosome().equals(chromosome()) && other.end() > start() && other.start() < end();
    }

    @Override
    default ChrBaseRegion chrBaseRegion()
    {
        return new ChrBaseRegion(chromosome(), start(), end());
    }
}
