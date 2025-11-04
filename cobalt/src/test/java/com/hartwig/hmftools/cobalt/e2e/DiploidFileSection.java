package com.hartwig.hmftools.cobalt.e2e;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public record DiploidFileSection(HumanChromosome chromosome, int start, int stop) implements ChromosomeSection
{
    @Override
    public String line(final int position)
    {
        return format("chr%s\t%d\t%d", chromosome().contig(),start,stop);
    }

    @Override
    public boolean is1Based()
    {
        return true;
    }

    @Override
    public List<String> lines()
    {
        return List.of(line(0));
    }
}
