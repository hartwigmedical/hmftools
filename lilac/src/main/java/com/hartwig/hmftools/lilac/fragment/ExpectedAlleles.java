package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.OptionalInt;

import org.apache.commons.lang3.NotImplementedException;

public class ExpectedAlleles
{
    private final int mDefaultAlleleCount;
    private final int[] mAlleleCount_;

    public ExpectedAlleles(final int defaultAlleleCount, final int[] alleleCount_)
    {
        mDefaultAlleleCount = defaultAlleleCount;
        mAlleleCount_ = alleleCount_;
    }

    private int expectedAlleles(int loci)
    {
        return mAlleleCount_ == null || loci >= mAlleleCount_.length ? mDefaultAlleleCount : mAlleleCount_[loci];
    }

    public int expectedAlleles(final List<Integer> loci)
    {
        return loci.stream().mapToInt(x -> expectedAlleles(x)).min().orElse(0);
    }

    public static ExpectedAlleles expectedAlleles(int geneCount, final NavigableMap<Integer, Integer> minUniqueProteinExonBoundaries_)
    {
        if(minUniqueProteinExonBoundaries_.isEmpty())
            return new ExpectedAlleles(2 * geneCount, null);

        int max = minUniqueProteinExonBoundaries_.lastKey();
        int arraySize = max - 1;
        int[] alleleCounts = new int[arraySize];

        int boundaryCount;
        for(int i = 0; i < arraySize; i++)
        {
            int locus = i + 1;
            boundaryCount = minUniqueProteinExonBoundaries_.headMap(locus, true).values().stream().mapToInt(x -> x).sum();
            alleleCounts[i] = 2 * (geneCount - boundaryCount);
        }

        boundaryCount = minUniqueProteinExonBoundaries_.values().stream().mapToInt(x -> x).sum();
        int defaultAlleleCount = 2 * (geneCount - boundaryCount);
        return new ExpectedAlleles(defaultAlleleCount, alleleCounts);
    }
}
