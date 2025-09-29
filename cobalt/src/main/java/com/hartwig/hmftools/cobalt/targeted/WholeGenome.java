package com.hartwig.hmftools.cobalt.targeted;

import com.hartwig.hmftools.cobalt.calculations.DoNothingNormaliser;
import com.hartwig.hmftools.cobalt.calculations.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class WholeGenome implements CobaltScope
{
    @Override
    public ResultsNormaliser finalNormaliser()
    {
        return new DoNothingNormaliser();
    }

    @Override
    public double enrichmentQuotient(final Chromosome chromosome, final DepthReading readDepth)
    {
        return 1.0;
    }

    @Override
    public boolean onTarget(final Chromosome chromosome, final int position)
    {
        return true; //todo
    }
}
