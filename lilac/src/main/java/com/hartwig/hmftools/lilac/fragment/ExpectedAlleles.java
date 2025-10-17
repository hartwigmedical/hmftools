package com.hartwig.hmftools.lilac.fragment;

import java.util.Collection;
import java.util.NavigableMap;

public class ExpectedAlleles
{
    private final int mDefaultAlleleCount;
    private final int[] mAlleleCount;

    public ExpectedAlleles(int defaultAlleleCount, final int[] alleleCount)
    {
        mDefaultAlleleCount = defaultAlleleCount;
        mAlleleCount = alleleCount;
    }

    private int expectedAlleles(int loci)
    {
        return mAlleleCount == null || loci >= mAlleleCount.length ? mDefaultAlleleCount : mAlleleCount[loci];
    }

    public int expectedAlleles(final Collection<Integer> loci)
    {
        return loci.stream().mapToInt(this::expectedAlleles).min().orElse(0);
    }

    public static ExpectedAlleles expectedAlleles(int geneCount, final NavigableMap<Integer, Integer> minUniqueProteinExonBoundaries)
    {
        if(minUniqueProteinExonBoundaries.isEmpty())
            return new ExpectedAlleles(2 * geneCount, null);

        int max = minUniqueProteinExonBoundaries.lastKey();
        int arraySize = max - 1;
        int[] alleleCounts = new int[arraySize];

        int boundaryCount;
        for(int i = 0; i < arraySize; i++)
        {
            int locus = i + 1;
            boundaryCount = minUniqueProteinExonBoundaries.headMap(locus, true).values().stream().mapToInt(x -> x).sum();
            alleleCounts[i] = 2 * (geneCount - boundaryCount);
        }

        boundaryCount = minUniqueProteinExonBoundaries.values().stream().mapToInt(x -> x).sum();
        int defaultAlleleCount = 2 * (geneCount - boundaryCount);
        return new ExpectedAlleles(defaultAlleleCount, alleleCounts);
    }
}
