package com.hartwig.hmftools.common.cobalt;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public interface CobaltDepth extends GenomePosition
{
    double referenceReadDepth();

    double tumorReadDepth();
}
