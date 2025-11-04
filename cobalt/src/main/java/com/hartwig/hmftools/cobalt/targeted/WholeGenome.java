package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import com.hartwig.hmftools.cobalt.calculations.DoNothingNormaliser;
import com.hartwig.hmftools.cobalt.calculations.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.calculations.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.consolidation.LowCoverageConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class WholeGenome implements CobaltScope
{
    @Override
    public ResultsNormaliser finalNormaliser()
    {
        return new DoNothingNormaliser();
    }

    @Override
    public ReadDepthStatisticsNormaliser medianByMeanNormaliser()
    {
        return new ReadDepthStatisticsNormaliser();
    }

    @Override
    public double enrichmentQuotient(final Chromosome chromosome, final int position)
    {
        return 1.0;
    }

    @Override
    public boolean onTarget(final Chromosome chromosome, final int position)
    {
        return true; //todo
    }

    @Override
    public ResultsConsolidator resultsConsolidator(final double medianReadDepth)
    {
        if(Double.isNaN(medianReadDepth))
        {
            return new NoOpConsolidator();
        }
        int consolidationCount = ResultsConsolidator.calcConsolidationCount(medianReadDepth);
        if(consolidationCount == 1)
        {
            CB_LOGGER.info("median read depth: {}, not using sparse consolidation", medianReadDepth);
            return new NoOpConsolidator();
        }
        CB_LOGGER.info("median read depth: {}, sparse consolidation count: {}", medianReadDepth, consolidationCount);
        return new LowCoverageConsolidator(consolidationCount);
    }
}
