package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

public class ExpectedAlleles
{
    private final int[] mAlleleCount;

    public ExpectedAlleles(final int[] alleleCount)
    {
        mAlleleCount = alleleCount;
    }

    public final int expectedAlleles(int loci)
    {
        return loci >= mAlleleCount.length ? 2 : mAlleleCount[loci];
    }

    public final int expectedAlleles(final List<Integer> loci)
    {
        // return loci.map { expectedAlleles(it) }.min()!!
        // CHECK:
        return loci.stream().mapToInt(x -> expectedAlleles(x)).min().orElse(0);
    }

    public static final ExpectedAlleles expectedAlleles(int otherMin1, int otherMin2)
    {
        /*
            fun expectedAlleles(otherMin1: Int, otherMin2: Int): ExpectedAlleles {
            val min = min(otherMin1, otherMin2)
            val max = max(otherMin1, otherMin2)

            val three = (1 until min).map { 6 }
            val two = (min until max).map { 4 }
            return ExpectedAlleles((three + two).toIntArray())
        }
         */

        // CHECK: seems to be 1 short on the 6-groups

        int min = min(otherMin1, otherMin2);
        int max = max(otherMin1, otherMin2);

        int arraySize = max - 1;
        int[] alleleCounts = new int[arraySize];

        for(int i = 0; i < min - 1; ++i)
        {
            alleleCounts[i] = 6;
        }

        for(int i = min - 1; i < max - 1; ++i)
        {
            alleleCounts[i] = 4;
        }

        return new ExpectedAlleles(alleleCounts);
    }
}
