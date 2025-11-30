package com.hartwig.hmftools.common.segmentation;

import java.util.List;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.segmentation.copynumber.PiecewiseConstantFit;

class AbsoluteSegments<T extends GenomePosition> extends ChromosomeArmSegments<T>
{
    AbsoluteSegments(ChrArm arm, final List<T> ratios, double[] values, PiecewiseConstantFit segmentation)
    {
        super(arm, ratios, values, segmentation);
    }

    public int segmentEnd(T endRatio)
    {
        return endRatio.position();
    }
}
