package com.hartwig.hmftools.cobalt.segmentation;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;

class PCFCoordinates
{
    private final PiecewiseConstantFit PCF;
    private final int WindowSize;
    private final int PositionOfFirstInterval;
    private final String Chromosome;

    PCFCoordinates(final PiecewiseConstantFit pcf, final int windowSize, int positionOfFirstInterval, final String chromosome)
    {
        PCF = pcf;
        WindowSize = windowSize;
        PositionOfFirstInterval = positionOfFirstInterval;
        Chromosome = chromosome;
    }

    List<ChrBaseRegion> intervals()
    {
        List<ChrBaseRegion> result = new ArrayList<>(PCF.lengths().length);
        for(int i = 0; i < PCF.lengths().length; i++)
        {
            int startPosition = PositionOfFirstInterval + PCF.startPositions()[i] * WindowSize;
            int endPosition = startPosition + PCF.lengths()[i] * WindowSize - 1;
            result.add(new ChrBaseRegion(Chromosome, startPosition, endPosition));
        }
        return result;
    }
}
