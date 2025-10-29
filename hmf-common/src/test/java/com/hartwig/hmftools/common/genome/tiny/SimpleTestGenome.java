package com.hartwig.hmftools.common.genome.tiny;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public abstract class SimpleTestGenome implements RefGenomeInterface
{
    @Override
    public String getBaseString(final String chromosome, final List<int[]> baseRanges)
    {
        StringBuilder refBases = new StringBuilder();
        baseRanges.forEach(x -> refBases.append(getBaseString(chromosome, x[0], x[1])));
        return refBases.toString();
    }

    @Override
    public int getChromosomeLength(final String chromosome)
    {
        return 0;
    }

    @Override
    public Map<String, Integer> chromosomeLengths()
    {
        return Map.of();
    }
}
