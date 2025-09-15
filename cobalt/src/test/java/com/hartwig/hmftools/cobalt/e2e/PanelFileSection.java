package com.hartwig.hmftools.cobalt.e2e;

import static java.lang.String.*;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public record PanelFileSection(HumanChromosome chromosome, int start, int stop, double relativeEnrichment) implements ChromosomeSection
{
    @Override
    public String line(final int position)
    {
        return format("chr%s\t%d\t%.9f", chromosome().toString(),position, relativeEnrichment);
    }

    @Override
    public boolean is1Based()
    {
        return true;
    }
}
