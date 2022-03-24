package com.hartwig.hmftools.common.cobalt;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public interface CobaltCount extends GenomePosition
{
    int referenceReadCount();

    int tumorReadCount();
}
