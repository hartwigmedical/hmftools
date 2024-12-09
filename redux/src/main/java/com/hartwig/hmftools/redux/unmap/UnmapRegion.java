package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class UnmapRegion extends ChrBaseRegion
{
    public final boolean IsStandardChromosome;
    public final int MaxDepth;

    public UnmapRegion(final String chromosome, final int posStart, final int posEnd, final int maxDepth)
    {
        super(chromosome, posStart, posEnd);
        IsStandardChromosome = HumanChromosome.contains(chromosome);
        MaxDepth = maxDepth;
    }

    public String toString()
    {
        return format("%s depth(%d)", super.toString(), MaxDepth);
    }
}
