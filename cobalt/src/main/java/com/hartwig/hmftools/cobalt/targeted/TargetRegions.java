package com.hartwig.hmftools.cobalt.targeted;

import com.hartwig.hmftools.cobalt.calculations.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public interface TargetRegions
{

    double enrichmentQuotient(Chromosome chromosome, DepthReading readDepth);

    boolean onTarget(Chromosome chromosome, int position);

    ResultsNormaliser createNormaliser();
}
