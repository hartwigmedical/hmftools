package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class SigUtils
{
    public static final Logger SU_LOGGER = LogManager.getLogger(SigUtils.class);

    public static void convertToPercentages(double[] counts)
    {
        double total = sumVector(counts);

        if(total <= 0)
            return;

        for(int i = 0; i < counts.length; ++i)
        {
            counts[i] /= total;
        }
    }

    public static double[] calculateFittedCounts(final Matrix signatures, final double[] allocations)
    {
        double[] fittedCounts = new double[signatures.Rows];

        for(int transId = 0; transId < signatures.Cols; ++transId)
        {
            double allocation = allocations[transId];

            for(int catId = 0; catId < signatures.Rows; ++catId)
            {
                fittedCounts[catId] += allocation * signatures.get(catId, transId);
            }
        }

        return fittedCounts;
    }

    public static SigResiduals calcResiduals(final double[] counts, final double[] fittedCounts, double totalCounts)
    {
        SigResiduals residuals = new SigResiduals();

        for(int catId = 0; catId < counts.length; ++catId)
        {
            double diff = fittedCounts[catId] - counts[catId];
            residuals.Total += abs(diff);

            if(diff > 0)
                residuals.Excess += diff;
        }

        residuals.Percent = residuals.Total / totalCounts;
        return residuals;
    }

    public static double calcLinearLeastSquares(final double[] params, final double[] data)
    {
        if(data.length != params.length)
            return 0;

        // returns the best ratio applying the params to the data assuming direct ratio (ie line through origin)
        double paramTotal = 0;
        double multTotal = 0;

        for(int i = 0; i < data.length; ++i)
        {
            paramTotal += params[i] * params[i];
            multTotal += params[i] * data[i];
        }

        return paramTotal > 0 ? multTotal/paramTotal : 0;
    }



    public static double calcAbsDiffs(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length)
            return 0;

        double diffTotal = 0;
        for(int i = 0; i < set1.length; ++i)
        {
            diffTotal += abs(set1[i] - set2[i]);
        }

        return diffTotal;
    }

}
