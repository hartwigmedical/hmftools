package com.hartwig.hmftools.cobalt.e2e;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

record CobaltRatioLine(String chromosome, int position, double tumorGcRatio) implements GenomePosition
{

}
