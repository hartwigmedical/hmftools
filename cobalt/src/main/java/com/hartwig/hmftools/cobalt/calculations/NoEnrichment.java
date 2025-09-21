package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class NoEnrichment implements CobaltCalculation.TargetRegions
{
    @Override
    public double enrichmentQuotient(final Chromosome chromosome, final ReadDepth readDepth)
    {
        return 1.0;
    }

    @Override
    public boolean isInTargetRegions(final Chromosome chromosome, final CobaltWindow window)
    {
        return true;
    }

    @Override
    public boolean isInTargetRegions(final Chromosome chromosome, final int position)
    {
        return true; //todo
    }

    @Override
    public boolean applyFinalNormalisation()
    {
        return false;
    }
}
