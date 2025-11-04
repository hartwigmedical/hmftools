package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public interface GenomeFilter
{
    boolean exclude(final Chromosome chromosome, DepthReading readDepth);
    Double referenceGcValueForWindow(final Chromosome chromosome, int position);
}
