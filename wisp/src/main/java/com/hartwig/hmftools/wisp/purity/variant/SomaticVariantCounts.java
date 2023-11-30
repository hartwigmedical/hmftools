package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

public class SomaticVariantCounts
{
    private int mAlleleFragments;
    private int mTotalFragments;
    private double mAllelelQualTotal;

    public SomaticVariantCounts()
    {
        mAlleleFragments = 0;
        mTotalFragments = 0;
        mAllelelQualTotal = 0;
    }

    public void addFragmentCount(int count)
    {
        mTotalFragments += count;
    }

    public void addAlleleFragmentCount(int count, double qualTotal)
    {
        mAlleleFragments += count;
        mAllelelQualTotal += qualTotal;
    }

    public int alleleFragments() { return mAlleleFragments; }
    public int totalFragments() { return mTotalFragments; }
    public double allelelQualTotal() { return mAllelelQualTotal; }

    public String toString()
    {
        return format("AF(%d) avgQualTotal(%.1f)", mAlleleFragments, mAllelelQualTotal);
    }
}
