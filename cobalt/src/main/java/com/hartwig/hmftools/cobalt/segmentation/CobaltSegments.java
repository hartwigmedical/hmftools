package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.common.utils.Doubles.mean;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;
import com.hartwig.hmftools.common.utils.pcf.CobaltSegment;

class CobaltSegments
{
    public final List<CobaltSegment> Segments = new ArrayList<>();
    public final ChrArm mArm;

    CobaltSegments(ChrArm arm, List<CobaltRatio> ratios, double[] values, PiecewiseConstantFit segmentation)
    {
        mArm = arm;
        String chr = mArm.chromosome().shortName();
        int ratiosIndex = 0;
        for(int i = 0; i < segmentation.lengths().length; i++)
        {
            int numberOfRatiosInSegment = segmentation.lengths()[i];
            double[] segmentValues = new double[numberOfRatiosInSegment];
            System.arraycopy(values, ratiosIndex, segmentValues, 0, numberOfRatiosInSegment);
            CobaltRatio startRatio = ratios.get(ratiosIndex);
            ratiosIndex = ratiosIndex + numberOfRatiosInSegment;
            CobaltRatio endRatio = ratios.get(ratiosIndex - 1);
            Segments.add(new CobaltSegment(chr, startRatio.position(), endRatio.position() + WINDOW_SIZE - 1, mean(segmentValues)));
        }
    }
}
