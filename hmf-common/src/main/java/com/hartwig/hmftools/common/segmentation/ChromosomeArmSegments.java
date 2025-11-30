package com.hartwig.hmftools.common.segmentation;

import static com.hartwig.hmftools.common.utils.Doubles.mean;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.segmentation.copynumber.PiecewiseConstantFit;
import com.hartwig.hmftools.common.utils.pcf.PcfSegment;

public abstract class ChromosomeArmSegments<T extends GenomePosition>
{
    public final List<PcfSegment> Segments = new ArrayList<>();
    public final ChrArm mArm;

    ChromosomeArmSegments(ChrArm arm, List<T> ratios, double[] values, PiecewiseConstantFit segmentation)
    {
        mArm = arm;
        String chr = mArm.chromosome().shortName();
        int ratiosIndex = 0;
        for(int i = 0; i < segmentation.lengths().length; i++)
        {
            int numberOfRatiosInSegment = segmentation.lengths()[i];
            double[] segmentValues = new double[numberOfRatiosInSegment];
            System.arraycopy(values, ratiosIndex, segmentValues, 0, numberOfRatiosInSegment);
            T startRatio = ratios.get(ratiosIndex);
            ratiosIndex = ratiosIndex + numberOfRatiosInSegment;
            T endRatio = ratios.get(ratiosIndex - 1);
            Segments.add(new PcfSegment(chr, startRatio.position(), segmentEnd(endRatio), mean(segmentValues)));
        }
    }

    public abstract int segmentEnd(T endRatio);
}

