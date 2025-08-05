package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

import org.apache.commons.lang3.NotImplementedException;

public class ExpectedAlleles
{
    private final int[] mAlleleCount_;

    public ExpectedAlleles(final int[] alleleCount_)
    {
        mAlleleCount_ = alleleCount_;
    }

    private int expectedAlleles(int loci)
    {
        return mAlleleCount_ == null || loci >= mAlleleCount_.length ? 2 : mAlleleCount_[loci];
    }

    public int expectedAlleles(final List<Integer> loci)
    {
        return loci.stream().mapToInt(x -> expectedAlleles(x)).min().orElse(0);
    }

    public static final ExpectedAlleles expectedAlleles(int otherMin1, int otherMin2)
    {
        // seems to be 1 short on the 6-groups
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
