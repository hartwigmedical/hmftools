package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.pow;

public class FragmentTotals
{
    private int mVariantCount;
    private int mTumorAdTotal;
    private int mSampleAdTotal;
    private double mSampleAllelelQualTotal;
    private double mTumorCopyNumberTotal;

    private int mTumorDepthTotal;
    private int mSampleDepthTotal;

    private int mSampleOneFragmentCount; // count of variants with 1 observed fragment
    private int mSampleTwoPlusCount;

    private double mSampleAdjustedAdTotal;
    private double mSampleAdjustedDepthTotal;

    private double mSampleWeightedDepthTotal;
    private double mSampleTumorAdjustedDepthTotal;

    private double mTumorAdjustedAdTotal;
    private double mTumorAdjustedDepthTotal;

    private Double mTumorVafOverride; // for ref variant analysis without a tumor sample

    public FragmentTotals()
    {
        mVariantCount = 0;
        mTumorAdTotal = 0;
        mSampleAdTotal = 0;
        mSampleAllelelQualTotal = 0;
        mTumorCopyNumberTotal = 0;
        mTumorDepthTotal = 0;
        mSampleDepthTotal = 0;
        mSampleOneFragmentCount = 0;
        mSampleTwoPlusCount = 0;
        mSampleAdjustedAdTotal = 0;
        mSampleAdjustedDepthTotal = 0;
        mSampleWeightedDepthTotal = 0;
        mSampleTumorAdjustedDepthTotal = 0;
        mTumorAdjustedAdTotal = 0;
        mTumorAdjustedDepthTotal = 0;
        mTumorVafOverride = null;
    }

    public void addVariantData(
            double copyNumber, int tumorAlleleFrags, int sampleAlleleFrags, int tumorDepth, int sampleDepth, double sampleQualTotal)
    {
        ++mVariantCount;
        mTumorAdTotal += tumorAlleleFrags;
        mTumorDepthTotal += tumorDepth;

        mSampleAdTotal += sampleAlleleFrags;
        mSampleDepthTotal += sampleDepth;
        mSampleAllelelQualTotal += sampleQualTotal;

        if(sampleAlleleFrags >= 2)
            ++mSampleTwoPlusCount;
        else if(sampleAlleleFrags == 1)
            ++mSampleOneFragmentCount;

        mTumorCopyNumberTotal += copyNumber;

        // wVAF = Σ(i=1->n)[ADi /CNi] * Σ(i=1->n)[DPi /CNi]

        mTumorAdjustedAdTotal += tumorAlleleFrags / copyNumber;
        mSampleAdjustedAdTotal += sampleAlleleFrags / copyNumber;

        double tumorDpPerCn = tumorDepth / copyNumber;
        mTumorAdjustedDepthTotal += tumorDpPerCn;

        mSampleAdjustedDepthTotal += sampleDepth / copyNumber;
        mSampleTumorAdjustedDepthTotal += sampleDepth * tumorDpPerCn;

        mSampleWeightedDepthTotal += tumorDpPerCn * pow(sampleDepth, 2);
    }

    public int variantCount() { return mVariantCount; }
    public int tumorAdTotal() { return mTumorAdTotal; }
    public int sampleAdTotal() { return mSampleAdTotal; }
    public int tumorDepthTotal() { return mTumorDepthTotal; }
    public int sampleDepthTotal() { return mSampleDepthTotal; }
    public int sampleOneFragmentCount() { return mSampleOneFragmentCount; }
    public int sampleTwoPlusCount() { return mSampleTwoPlusCount; }

    public void setTumorVafOverride(final double vaf) { mTumorVafOverride = vaf; }

    public double rawTumorVaf()
    {
        if(mTumorVafOverride != null)
            return mTumorVafOverride;

        return mTumorDepthTotal > 0 ? mTumorAdTotal / (double)mTumorDepthTotal : 0;
    }

    public double adjTumorVaf()
    {
        if(mTumorVafOverride != null)
            return mTumorVafOverride;

        return mTumorAdjustedDepthTotal > 0 ? mTumorAdjustedAdTotal / mTumorAdjustedDepthTotal : 0;
    }

    public double adjSampleVaf() { return mSampleAdjustedDepthTotal > 0 ? mSampleAdjustedAdTotal / mSampleAdjustedDepthTotal : 0; }

    public double adjSampleVaf(double sampleAdAdjustment)
    {
        if(mSampleAdjustedDepthTotal == 0)
            return 0;

        double avgCopyNumber = mTumorCopyNumberTotal / mVariantCount;
        double adjSampleAdTotal = mSampleAdjustedAdTotal + sampleAdAdjustment / avgCopyNumber;

        return adjSampleAdTotal / mSampleAdjustedDepthTotal;
    }

    public double weightedSampleDepth()
    {
        // wAD = Σ(i=1->n)[(DPi_cfDNA)^2*DPi _tissue/CNn] / Σ(i=1->n)[DPi_cfDNA *DPi_Tissue/CNi]
        return mSampleTumorAdjustedDepthTotal > 0 ? mSampleWeightedDepthTotal / mSampleTumorAdjustedDepthTotal : 0;
    }

    public double qualPerAlleleFragment()
    {
        return mSampleAdTotal > 0 ? mSampleAllelelQualTotal / mSampleAdTotal : 0;
    }
}
