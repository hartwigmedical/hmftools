package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.pow;

public class FragmentTotals
{
    private int mVariantCount;
    private int mTumorFragsTotal;
    private int mSampleFragsTotal;
    private double mSampleAllelelQualTotal;
    private double mTumorCopyNumberTotal;

    private int mTumorDepthTotal;
    private int mSampleDepthTotal;

    private int mSampleOneFragmentCount; // count of variants with 1 observed fragment
    private int mSampleTwoPlusCount;

    private double mSampleAdjustedFragsTotal;
    private double mSampleAdjustedDepthTotal;

    private double mSampleWeightedDepthTotal;
    private double mSampleTumorAdjustedDepthTotal;
    private double mSampleWeightedFragsTotal;

    private double mTumorAdjustedFragsTotal;
    private double mTumorAdjustedDepthTotal;

    private Double mTumorVafOverride; // for ref variant analysis without a tumor sample

    public FragmentTotals()
    {
        mVariantCount = 0;
        mTumorFragsTotal = 0;
        mSampleFragsTotal = 0;
        mSampleAllelelQualTotal = 0;
        mTumorCopyNumberTotal = 0;
        mTumorDepthTotal = 0;
        mSampleDepthTotal = 0;
        mSampleOneFragmentCount = 0;
        mSampleTwoPlusCount = 0;
        mSampleAdjustedFragsTotal = 0;
        mSampleAdjustedDepthTotal = 0;
        mSampleWeightedDepthTotal = 0;
        mSampleTumorAdjustedDepthTotal = 0;
        mSampleWeightedFragsTotal = 0;
        mTumorAdjustedFragsTotal = 0;
        mTumorAdjustedDepthTotal = 0;
        mTumorVafOverride = null;
    }

    public void addVariantData(
            double copyNumber, int tumorAlleleFrags, int sampleAlleleFrags, int tumorDepth, int sampleDepth, double sampleQualTotal)
    {
        ++mVariantCount;
        mTumorFragsTotal += tumorAlleleFrags;
        mTumorDepthTotal += tumorDepth;

        mSampleFragsTotal += sampleAlleleFrags;
        mSampleDepthTotal += sampleDepth;
        mSampleAllelelQualTotal += sampleQualTotal;

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

        mSampleAdjustedDepthTotal += sampleDepth / copyNumber;
        mSampleTumorAdjustedDepthTotal += sampleDepth * tumorDpPerCn;

        mSampleWeightedFragsTotal += tumorDpPerCn * sampleAlleleFrags * sampleDepth;

        mSampleWeightedDepthTotal += tumorDpPerCn * pow(sampleDepth, 2);
    }

    public int variantCount() { return mVariantCount; }
    public int tumorAdTotal() { return mTumorFragsTotal; }
    public int sampleAdTotal() { return mSampleFragsTotal; }
    public int tumorDepthTotal() { return mTumorDepthTotal; }
    public int sampleDepthTotal() { return mSampleDepthTotal; }
    public double sampleAdjustedDepthTotal() { return mSampleAdjustedDepthTotal; }
    public double sampleTumorAdjustedDepthTotal() { return mSampleTumorAdjustedDepthTotal; }
    public double sampleTumorAdjustedFragsTotal() { return mSampleWeightedFragsTotal; }
    public int sampleOneFragmentCount() { return mSampleOneFragmentCount; }
    public int sampleTwoPlusCount() { return mSampleTwoPlusCount; }

    public void setTumorVafOverride(final double vaf) { mTumorVafOverride = vaf; }

    public double rawTumorVaf()
    {
        if(mTumorVafOverride != null)
            return mTumorVafOverride;

        return mTumorDepthTotal > 0 ? mTumorFragsTotal / (double)mTumorDepthTotal : 0;
    }

    public double adjTumorVaf()
    {
        if(mTumorVafOverride != null)
            return mTumorVafOverride;

        return mTumorAdjustedDepthTotal > 0 ? mTumorAdjustedFragsTotal / mTumorAdjustedDepthTotal : 0;
    }

    public double adjSampleVaf() { return mSampleAdjustedDepthTotal > 0 ? mSampleAdjustedFragsTotal / mSampleAdjustedDepthTotal : 0; }

    public double adjSampleVaf(double sampleAdAdjustment)
    {
        if(mSampleAdjustedDepthTotal == 0)
            return 0;

        double avgCopyNumber = mTumorCopyNumberTotal / mVariantCount;
        double adjSampleAdTotal = mSampleAdjustedFragsTotal + sampleAdAdjustment / avgCopyNumber;

        return adjSampleAdTotal / mSampleAdjustedDepthTotal;
    }

    public double weightedSampleDepth()
    {
        // wAD = Σ(i=1->n)[(DPi_cfDNA)^2*DPi _tissue/CNn] / Σ(i=1->n)[DPi_cfDNA *DPi_Tissue/CNi]
        return mSampleTumorAdjustedDepthTotal > 0 ? mSampleWeightedDepthTotal / mSampleTumorAdjustedDepthTotal : 0;
    }

    public double weightedSampleFrags()
    {
        // wAD = Σ(i=1->n)[(DPi_cfDNA)^2*DPi _tissue/CNn] / Σ(i=1->n)[DPi_cfDNA *DPi_Tissue/CNi]
        return mSampleTumorAdjustedDepthTotal > 0 ? mSampleWeightedFragsTotal / mSampleTumorAdjustedDepthTotal : 0;
    }

    public double qualPerAlleleFragment()
    {
        return mSampleFragsTotal > 0 ? mSampleAllelelQualTotal / mSampleFragsTotal : 0;
    }
}
