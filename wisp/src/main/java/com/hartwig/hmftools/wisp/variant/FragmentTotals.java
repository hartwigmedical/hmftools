package com.hartwig.hmftools.wisp.variant;

import static com.hartwig.hmftools.wisp.WispConstants.MAX_COPY_NUMBER_FOR_FITTING;
import static java.lang.Math.max;

public class FragmentTotals
{
    private int mVariantCount;
    private int mTumorFragsTotal;
    private int mSampleFragsTotal;

    private int mTumorDepthTotal;
    private int mSampleDepthTotal;
    private double mDepthToCopyNumberWeightingTotal;
    private double mSampleWeightedAfTotal;

    private int mOneFragmentCount; // count of variants with 1 observed fragment
    private int mTwoPlusCount;

    private double mWeightedOneFragmentCount;
    private double mWeightedTwoPlusCount;

    private double mSampleWeightedDepthTotal;

    private double mVcnSampleDepthTotal;
    private double mCnSampleDepthTotal;

    private Double mTumorVafOverride; // for ref variant analysis without a tumor sample

    public FragmentTotals()
    {
        mVariantCount = 0;
        mTumorFragsTotal = 0;
        mSampleFragsTotal = 0;
        mTumorDepthTotal = 0;
        mSampleDepthTotal = 0;
        mDepthToCopyNumberWeightingTotal = 0;
        mSampleWeightedAfTotal = 0;
        mOneFragmentCount = 0;
        mTwoPlusCount = 0;
        mWeightedOneFragmentCount = 0;
        mWeightedTwoPlusCount = 0;
        mSampleWeightedDepthTotal = 0;
        mVcnSampleDepthTotal = 0;
        mCnSampleDepthTotal = 0;

        mTumorVafOverride = null;
    }

    public void addVariantData(
            double copyNumber, double variantCopyNumber, int tumorAlleleFrags, int sampleAlleleFrags, int tumorDepth, int sampleDepth,
            boolean capMaxTumorCopyNumber)
    {
        if(capMaxTumorCopyNumber && copyNumber > MAX_COPY_NUMBER_FOR_FITTING)
            return;

        ++mVariantCount;
        mTumorFragsTotal += tumorAlleleFrags;
        mTumorDepthTotal += tumorDepth;

        mSampleFragsTotal += sampleAlleleFrags;
        mSampleDepthTotal += sampleDepth;

        double depthToCopyNumberFactor = sampleDepth / max(copyNumber, 1);  // ratio of actual to expected sample DP (low for off-target variants)
        mDepthToCopyNumberWeightingTotal += depthToCopyNumberFactor;

        if(sampleAlleleFrags >= 2)
        {
            ++mTwoPlusCount;
            mWeightedTwoPlusCount += depthToCopyNumberFactor;
        }
        else if(sampleAlleleFrags == 1)
        {
            ++mOneFragmentCount;
            mWeightedOneFragmentCount += depthToCopyNumberFactor;
        }

        if(sampleDepth > 0)
            mSampleWeightedAfTotal += (sampleAlleleFrags / (double)sampleDepth) * depthToCopyNumberFactor;

        mSampleWeightedDepthTotal += sampleDepth * depthToCopyNumberFactor; // numerator for WAD

        // for WA_VCN and WA_CN
        mVcnSampleDepthTotal += variantCopyNumber * depthToCopyNumberFactor;
        mCnSampleDepthTotal += copyNumber * depthToCopyNumberFactor;
    }

    public int variantCount() { return mVariantCount; }
    public int sampleAdTotal() { return mSampleFragsTotal; }
    public int sampleDepthTotal() { return mSampleDepthTotal; }
    public int oneFragmentCount() { return mOneFragmentCount; }
    public int twoPlusCount() { return mTwoPlusCount; }
    public double weightedOneFragmentCount() { return mWeightedOneFragmentCount; }
    public double weightedTwoPlusCount() { return mWeightedTwoPlusCount; }

    public void setTumorVafOverride(final double vaf) { mTumorVafOverride = vaf; }

    public double rawTumorVaf()
    {
        if(mTumorVafOverride != null)
            return mTumorVafOverride;

        return mTumorDepthTotal > 0 ? mTumorFragsTotal / (double)mTumorDepthTotal : 0;
    }

    public double rawSampleVaf() { return mSampleDepthTotal > 0 ? mSampleFragsTotal / (double)mSampleDepthTotal : 0; }

    public double adjSampleVaf() { return mDepthToCopyNumberWeightingTotal > 0 ? mSampleWeightedAfTotal / mDepthToCopyNumberWeightingTotal : 0; }

    public double weightedSampleDepth()
    {
        return mDepthToCopyNumberWeightingTotal > 0 ? mSampleWeightedDepthTotal / mDepthToCopyNumberWeightingTotal : 0;
    }

    // for WA_VCN and WA_CN
    public double weightedVariantCopyNumber()
    {
        return mDepthToCopyNumberWeightingTotal > 0 ? mVcnSampleDepthTotal / mDepthToCopyNumberWeightingTotal : 0;
    }

    public double weightedCopyNumber()
    {
        return mDepthToCopyNumberWeightingTotal > 0 ? mCnSampleDepthTotal / mDepthToCopyNumberWeightingTotal : 0;
    }
}
