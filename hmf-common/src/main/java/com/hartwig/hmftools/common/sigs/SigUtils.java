package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;

public class SigUtils
{
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

    public static final double[] calculateFittedCounts(final SigMatrix signatures, final double[] allocations)
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

    public static final int RESIDUAL_TOTAL = 0;
    public static final int RESIDUAL_PERC = 1;

    public static double[] calcResiduals(final double[] transcriptCounts, final double[] fittedCounts, double totalCounts)
    {
        double residualsTotal = 0;

        for(int catId = 0; catId < transcriptCounts.length; ++catId)
        {
            residualsTotal += abs(transcriptCounts[catId] - fittedCounts[catId]);
        }

        double residualsPerc = residualsTotal / totalCounts;

        return new double[] {residualsTotal, residualsPerc};

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

    public static void copyMatrix(final double[][] source, final double[][] dest)
    {
        if(source.length != dest.length)
            return;

        int rows = source.length;
        int cols = source[0].length;

        for(int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                dest[i][j] = source[i][j];
            }
        }
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
