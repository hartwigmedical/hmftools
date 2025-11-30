package com.hartwig.hmftools.common.segmentation;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.WINDOW_SIZE;

import java.util.List;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.segmentation.copynumber.PiecewiseConstantFit;

class WindowSegments<T extends GenomePosition> extends ChromosomeArmSegments<T>
{
    WindowSegments(ChrArm arm, final List<T> ratios, double[] values, PiecewiseConstantFit segmentation)
    {
        super(arm, ratios, values, segmentation);
    }

    public int segmentEnd(T endRatio)
    {
        return endRatio.position() + WINDOW_SIZE - 1;
    }
}
