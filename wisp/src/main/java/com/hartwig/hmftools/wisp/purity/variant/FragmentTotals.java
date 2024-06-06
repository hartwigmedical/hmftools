package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.pow;

public class FragmentTotals
{
    private int mVariantCount;
    private int mTumorFragsTotal;
    private int mSampleFragsTotal;
    private double mTumorCopyNumberTotal;

    private int mTumorDepthTotal;
    private int mSampleDepthTotal;

    private int mSampleOneFragmentCount; // count of variants with 1 observed fragment
    private int mSampleTwoPlusCount;

    private double mSampleAdjustedFragsTotal;

    private double mSampleDepthPerCopyNumberTotal;

    private double mSampleWeightedDepthTotal;
    private double mSampleTumorAdjustedDepthTotal;

    private double mTumorAdjustedFragsTotal;
    private double mTumorAdjustedDepthTotal;

    private double mVcnSampleDepthTotal;
    private double mCnSampleDepthTotal;

    private Double mTumorVafOverride; // for ref variant analysis without a tumor sample

    public FragmentTotals()
    {
        mVariantCount = 0;
        mTumorFragsTotal = 0;
        mSampleFragsTotal = 0;
        mTumorCopyNumberTotal = 0;
        mTumorDepthTotal = 0;
        mSampleDepthTotal = 0;
        mSampleOneFragmentCount = 0;
        mSampleTwoPlusCount = 0;
        mSampleAdjustedFragsTotal = 0;
        mSampleDepthPerCopyNumberTotal = 0;
        mSampleWeightedDepthTotal = 0;
        mSampleTumorAdjustedDepthTotal = 0;
        mTumorAdjustedFragsTotal = 0;
        mTumorAdjustedDepthTotal = 0;
        mVcnSampleDepthTotal = 0;
        mCnSampleDepthTotal = 0;

        mTumorVafOverride = null;
    }

    public void addVariantData(
            double copyNumber, double variantCopyNumber,
            int tumorAlleleFrags, int sampleAlleleFrags, int tumorDepth, int sampleDepth, double sampleQualTotal)
    {
        ++mVariantCount;
        mTumorFragsTotal += tumorAlleleFrags;
        mTumorDepthTotal += tumorDepth;

        mSampleFragsTotal += sampleAlleleFrags;
        mSampleDepthTotal += sampleDepth;

        if(sampleAlleleFrags >= 2)
            ++mSampleTwoPlusCount;
        else if(sampleAlleleFrags == 1)
            ++mSampleOneFragmentCount;

        mTumorCopyNumberTotal += copyNumber;

        // wVAF = Σ(i=1->n)[ADi /CNi] * Σ(i=1->n)[DPi /CNi]

        mTumorAdjustedFragsTotal += tumorAlleleFrags / copyNumber;
        mSampleAdjustedFragsTotal += sampleAlleleFrags / copyNumber;

        double tumorDpPerCn = tumorDepth / copyNumber;
        mTumorAdjustedDepthTotal += tumorDpPerCn;

        mSampleDepthPerCopyNumberTotal += sampleDepth / copyNumber;

        mSampleTumorAdjustedDepthTotal += sampleDepth * tumorDpPerCn; // denominator for weighted average depth (WAD)
        mSampleWeightedDepthTotal += tumorDpPerCn * pow(sampleDepth, 2); // numerator for WAD

        // for WA_VCN and WA_CN
        mVcnSampleDepthTotal += variantCopyNumber * sampleDepth;
        mCnSampleDepthTotal += copyNumber * sampleDepth;
    }

    public int variantCount() { return mVariantCount; }
    public int sampleAdTotal() { return mSampleFragsTotal; }
    public int tumorDepthTotal() { return mTumorDepthTotal; }
    public int sampleDepthTotal() { return mSampleDepthTotal; }
    public int sampleOneFragmentCount() { return mSampleOneFragmentCount; }
    public int sampleTwoPlusCount() { return mSampleTwoPlusCount; }

    public void setTumorVafOverride(final double vaf) { mTumorVafOverride = vaf; }

    public double rawTumorVaf()
    {
        if(mTumorVafOverride != null)
            return mTumorVafOverride;

        return mTumorDepthTotal > 0 ? mTumorFragsTotal / (double)mTumorDepthTotal : 0;
    }

    public double rawSampleVaf() { return mSampleDepthTotal > 0 ? mSampleFragsTotal / (double)mSampleDepthTotal : 0; }

    public double adjSampleVaf() { return mSampleDepthPerCopyNumberTotal > 0 ? mSampleAdjustedFragsTotal / mSampleDepthPerCopyNumberTotal : 0; }

    public double weightedSampleDepth()
    {
        return mSampleTumorAdjustedDepthTotal > 0 ? mSampleWeightedDepthTotal / mSampleTumorAdjustedDepthTotal : 0;
    }

    // for WA_VCN and WA_CN
    public double weightedVariantCopyNumber() { return mSampleDepthTotal > 0 ? mVcnSampleDepthTotal / mSampleDepthTotal : 0; }
    public double weightedCopyNumber() { return mSampleDepthTotal > 0 ? mCnSampleDepthTotal / mSampleDepthTotal : 0; }
}
