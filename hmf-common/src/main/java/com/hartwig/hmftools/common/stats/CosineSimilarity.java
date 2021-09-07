package com.hartwig.hmftools.common.stats;

import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

public final class CosineSimilarity
{
    public static double calcCosineSim(final double[] set1, final double[] set2)
    {
        return calcCosineSim(set1, set2, false, false);
    }

    public static double calcCosineSim(final double[] set1, final double[] set2, boolean skipZeros)
    {
        return calcCosineSim(set1, set2, skipZeros, skipZeros);
    }

    public static double calcCosineSim(final double[] set1, final double[] set2, boolean skipZeros1, boolean skipZeros2)
    {
        if(set1.length != set2.length || set1.length == 0)
            return 0;

        double aaTotal = 0;
        double bbTotal = 0;
        double abTotal = 0;
        int nonZeroCount = 0;

        for(int i = 0; i < set1.length; ++i)
        {
            double val1 = set1[i];
            double val2 = set2[i];

            if((skipZeros1 && val1 == 0) || (skipZeros2 && val2 == 0))
                continue;

            aaTotal += val1 * val1;
            bbTotal += val2 * val2;
            abTotal += val1 * val2;
            ++nonZeroCount;
        }

        if(aaTotal <= 0 || bbTotal <= 0 || nonZeroCount < 2)
            return 0;

        return min(abTotal / (sqrt(aaTotal) * sqrt(bbTotal)), 1.0);
    }

    public static double calcCSSRelative(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length || set1.length == 0)
            return 0;

        double aaTotal = 0;
        double bbTotal = 0;
        double abTotal = 0;

        // reduce each pair of values to be roughly in line with the median pair
        double total1 = sumVector(set1);
        double total2 = sumVector(set2);

        if(total1 == 0 || total2 == 0)
            return 0;

        double total = (total1 + total2) * 0.5;

        for(int i = 0; i < set1.length; ++i)
        {
            double a = set1[i];
            double b = set2[i];

            if(a == 0 || b == 0)
                continue;

            double ratioToTotal = ((a + b) * 0.5) / total;

            a /= ratioToTotal;
            b /= ratioToTotal;

            aaTotal += a*a;
            bbTotal += b*b;
            abTotal += a*b;
        }

        if(aaTotal <= 0 || bbTotal <= 0)
            return 0;

        return min(abTotal / (sqrt(aaTotal) * sqrt(bbTotal)), 1.0);
    }

}
