package com.hartwig.hmftools.cobalt.segmentation;

import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;

import org.apache.commons.math3.util.FastMath;

public class SegmentationData
{
    final private double[] valuesForSegmentation;
    final private double[] rawValues;

    public SegmentationData(List<CobaltRatio> ratios, Function<CobaltRatio, Double> value)
    {
        valuesForSegmentation = new double[ratios.size()];
        rawValues = new double[ratios.size()];
        for(int i = 0; i < ratios.size(); i++)
        {
            CobaltRatio ratio = ratios.get(i);
            final double v = value.apply(ratio);
            rawValues[i] = v;
            // Our R script that called copynumber put 0.001 as a floor for the ratios and converted
            // them to log_2 values. Note that negative ratio values have already been filtered out.
            if(v < 0.001)
            {
                valuesForSegmentation[i] = -9.965784;
            }
            else {
                valuesForSegmentation[i] = (float) FastMath.log(2, v);
            }
        }
    }

    public int count()
    {
        return valuesForSegmentation.length;
    }

    public boolean isEmpty()
    {
        return count() == 0;
    }

    public double[] valuesForSegmentation()
    {
        return valuesForSegmentation;
    }

    public double[] rawValues()
    {
        return rawValues;
    }
}
