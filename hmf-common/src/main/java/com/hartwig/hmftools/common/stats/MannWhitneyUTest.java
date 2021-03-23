package com.hartwig.hmftools.common.stats;

import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

public class MannWhitneyUTest
{
    private final NormalDistribution mStandardNormal;

    public MannWhitneyUTest()
    {
        mStandardNormal = new NormalDistribution((RandomGenerator)null, 0.0D, 1.0D);
    }

    public MwuResult calculate(final double[] values1, final double[] values2)
    {
        return calculate(values1, values2, mStandardNormal);
    }

    public static MwuResult calculate(final double[] values1, final double[] values2, final NormalDistribution normalDistribution)
    {
        int valueCount = values1.length + values2.length;
        final double[] combinedValues = new double[valueCount];
        final boolean[] isCohort1 = new boolean[valueCount];

        combineCohortValues(values1, values2, combinedValues, isCohort1);

        return calculate(combinedValues, isCohort1, normalDistribution);
    }

    public MwuResult calculate(final double[] values, final boolean[] isCohort1)
    {
        return calculate(values, isCohort1, mStandardNormal);
    }

    public static MwuResult calculate(final double[] values, final boolean[] isCohort1, final NormalDistribution normalDistribution)
    {
        if(values.length != isCohort1.length)
            return new MwuResult(false);

        MwuResult result = new MwuResult(true);

        // assign zero values equal ranks
        int zeroCount1 = 0;
        int zeroCount2 = 0;
        for(int i = 0; i < values.length; ++i)
        {
            if(values[i] > 0)
                break;

            if(isCohort1[i])
                ++zeroCount1;
            else
                ++zeroCount2;
        }

        if(zeroCount1 > 0 || zeroCount2 > 0)
        {
            int totalZeroValues = zeroCount1 + zeroCount2;
            int zeroRankSum = totalZeroValues * (totalZeroValues + 1) / 2;
            double cohortZeroValueRank1 = zeroCount1 / (double)totalZeroValues * zeroRankSum;
            result.RankSum1 += (int)round(cohortZeroValueRank1);
        }

        int count1 = zeroCount1;

        for(int i = 0; i < values.length; ++i)
        {
            if(values[i] == 0)
                continue;

            if(isCohort1[i])
            {
                result.RankSum1 += i + 1; // rank starts at 1 not 0
                ++count1;
            }
        }

        int count2 = values.length - count1;
        double n1n2 = count1 * count2;
        double mean = n1n2 * 0.5;
        double variance = sqrt(n1n2 * (count1 + count2 + 1) / 12.0);
        double u1 = result.RankSum1 - count1 * (count1 + 1) / 2.0;
        double u2 = n1n2 - u1;
        double uMin = min(u1, u2);
        double z = (uMin - mean) / variance;

        result.PValue = normalDistribution.cumulativeProbability(z) * 2;
        result.Count1 = count1;

        return result;
    }

    private static final double COHORT_1 = 1.0;
    private static final double COHORT_2 = 0.0;

    // take 2 distributions and form a single ascending set of values, with an accompanying array indicating the source of each
    public static void combineCohortValues(
            final double[] values1, final double[] values2, final double[] combinedValues, final boolean[] isCohort1)
    {
        if(combinedValues.length != values1.length + values2.length || isCohort1.length != combinedValues.length)
            return;

        // values in ascending order, with first array value being the value, second being the source (1/true for cohort 1)
        final List<double[]> combinedList = Lists.newArrayListWithExpectedSize(combinedValues.length);

        for(int i = 0; i <= 1; ++i)
        {
            final double[] values = (i == 0) ? values1 : values2;
            double isC1 = (i == 0) ? COHORT_1 : COHORT_2;

            for(int j = 0; j < values.length; ++j)
            {
                int index = 0;

                while(index < combinedList.size())
                {
                    if(values[j] < combinedList.get(index)[0])
                        break;
                    else
                        ++index;
                }

                combinedList.add(index, new double[] { values[j], isC1 });
            }
        }

        for(int i = 0; i < combinedValues.length; ++i)
        {
            final double[] pair = combinedList.get(i);
            combinedValues[i] = pair[0];
            isCohort1[i] = pair[1] == COHORT_1;
        }
    }
}
