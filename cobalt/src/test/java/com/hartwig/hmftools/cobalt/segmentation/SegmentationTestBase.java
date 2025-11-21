package com.hartwig.hmftools.cobalt.segmentation;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class SegmentationTestBase
{

    CobaltRatio ratioForValue(HumanChromosome chromosome, int start, double value)
    {
        // The RatioSegmenter uses the log of the values.
        // So that we know how the segmentation works out, we will anti-log them first.
        double v = value >= 0 ? Math.pow(2, value) : -1.0;
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }

    CobaltRatio ratio(HumanChromosome chromosome, int start, double v)
    {
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }
}
