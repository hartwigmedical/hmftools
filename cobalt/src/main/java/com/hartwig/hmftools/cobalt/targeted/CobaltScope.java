package com.hartwig.hmftools.cobalt.targeted;

import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public interface CobaltScope
{
    ResultsNormaliser finalNormaliser();

    ReadDepthStatisticsNormaliser medianByMeanNormaliser();

    ResultsConsolidator resultsConsolidator(double medianReadDepth);

    double findRegionEnrichment(final HumanChromosome chromosome, int position);

    boolean onTarget(final HumanChromosome chromosome, int position);
}
