package com.hartwig.hmftools.cobalt.targeted;

import com.hartwig.hmftools.cobalt.calculations.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.calculations.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public interface CobaltScope
{

    ResultsNormaliser finalNormaliser();

    ReadDepthStatisticsNormaliser medianByMeanNormaliser();

    ResultsConsolidator resultsConsolidator(double medianReadDepth);

    double enrichmentQuotient(Chromosome chromosome, int position);

    boolean onTarget(Chromosome chromosome, int position);
}
