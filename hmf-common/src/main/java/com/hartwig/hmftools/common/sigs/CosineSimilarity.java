package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.min;
import static java.lang.Math.sqrt;

public class CosineSimilarity
{
    public static double calcCosineSim(final double[] set1, final double[] set2)
    {
        return calcCosineSim(set1, set2, false);
    }

    public static double calcCosineSim(final double[] set1, final double[] set2, boolean skipZeros)
    {
        if(set1.length != set2.length || set1.length == 0)
            return 0;

        double aaTotal = 0;
        double bbTotal = 0;
        double abTotal = 0;
        int nonZeroCount = 0;

        for(int i = 0; i < set1.length; ++i)
        {
            double a = set1[i];
            double b = set2[i];

            if(skipZeros && (a == 0 || b == 0))
                continue;

            aaTotal += a*a;
            bbTotal += b*b;
            abTotal += a*b;
            ++nonZeroCount;
        }

        if(aaTotal <= 0 || bbTotal <= 0 || nonZeroCount < 2)
            return 0;

        return min(abTotal / (sqrt(aaTotal) * sqrt(bbTotal)), 1.0);
    }

}
