package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Integers.median;

import java.util.List;

import com.google.common.collect.Lists;

public class SomaticVariantCounts
{
    private int mAlleleFragments;
    private int mTotalFragments;
    private double mAllelelQualTotal;
    private final List<Integer> mVariantDepths;
    private final List<Integer> mNonZeroVariantDepths;

    public SomaticVariantCounts()
    {
        mVariantDepths = Lists.newArrayList();
        mNonZeroVariantDepths = Lists.newArrayList();
        mAlleleFragments = 0;
        mTotalFragments = 0;
        mAllelelQualTotal = 0;
    }

    public void addFragmentCount(int count)
    {
        mTotalFragments += count;
        mVariantDepths.add(count);

        if(count > 0)
            mNonZeroVariantDepths.add(count);
    }

    public void addAlleleFragmentCount(int count, double qualTotal)
    {
        mAlleleFragments += count;
        mAllelelQualTotal += qualTotal;
    }

    public int alleleFragments() { return mAlleleFragments; }
    public int totalFragments() { return mTotalFragments; }
    public double allelelQualTotal() { return mAllelelQualTotal; }

    public int nonZeroDepth() { return mNonZeroVariantDepths.size(); }

    public double medianDepth(boolean useNonZero)
    {
        return useNonZero ? median(mNonZeroVariantDepths) : median(mVariantDepths);
    }

    public String toString()
    {
        return format("AF(%d) avgQualTotal(%.1f) depthCounts(nonzero=%d all=%d)",
                mAlleleFragments, mAllelelQualTotal, mNonZeroVariantDepths.size(), mVariantDepths.size());
    }
}
