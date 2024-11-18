package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import com.hartwig.hmftools.wisp.purity.PurityConstants;

public class GenotypeFragments
{
    public final String SampleName;
    public final int AlleleCount;
    public final int Depth;
    public final double QualTotal;

    // where UMI counts are available then the Depth and AlleleCount are just the total from these, otherwise they are standard DP and AD
    public final UmiTypeCounts UmiCounts;

    private double mBqrErrorRate; // only used for output data

    private boolean mLikelyChip;

    public GenotypeFragments(final String sampleName, final int alleleCount, final int depth, final double qualTotal,
            final UmiTypeCounts umiCounts)
    {
        SampleName = sampleName;
        AlleleCount = alleleCount;
        Depth = depth;
        QualTotal = qualTotal;
        UmiCounts = umiCounts;

        mBqrErrorRate = 0;
        mLikelyChip = false;
    }

    public double qualPerAlleleFragment() { return AlleleCount > 0 ? QualTotal / (double)AlleleCount : 0; }

    public double vaf() { return Depth > 0 ? AlleleCount / (double)Depth : 0; }

    public boolean isLowQual() { return AlleleCount > 0 && qualPerAlleleFragment() <= PurityConstants.MIN_QUAL_PER_AD; }  // also require that MQF > 10

    public void setBqrErrorRate(double errorRate) { mBqrErrorRate = errorRate; }
    public double bqrErrorRate() { return mBqrErrorRate; }

    public void markLikeChip() { mLikelyChip = true; }
    public boolean likelyChip() { return mLikelyChip; }

    public String toString()
    {
        return format("%d/%d umi(total=%d allele=%d)", AlleleCount, Depth, UmiCounts.totalCount(), UmiCounts.alleleCount());
    }
}
