package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;

import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

import htsjdk.variant.variantcontext.Genotype;

public class GenotypeFragments
{
    public final String SampleName;
    public final int AlleleCount;
    public final int Depth;
    public final double QualTotal;
    public final Genotype GenotypeData;

    // where UMI counts are available then the Depth and AlleleCount are just the total from these, otherwise they are standard DP and AD
    public final UmiTypeCounts UmiCounts;

    private double mBqrErrorRate; // only used for output data

    private boolean mLikelyChip;

    public GenotypeFragments(
            final String sampleName, final int alleleCount, final int depth, final double qualTotal, final UmiTypeCounts umiCounts,
            final Genotype genotype)
    {
        SampleName = sampleName;
        AlleleCount = alleleCount;
        Depth = depth;
        QualTotal = qualTotal;
        UmiCounts = umiCounts;
        GenotypeData = genotype;

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

    public int averageReadDistance()
    {
        Object intField = GenotypeData.getExtendedAttribute(SageVcfTags.AVG_READ_EDGE_DISTANCE, null);

        if(intField == null)
            return -1;

        String[] aedCounts = intField.toString().split(LIST_SEPARATOR, 2);
        return aedCounts.length == 2 ? Integer.parseInt(aedCounts[1]) : -1;
    }

    public String toString()
    {
        return format("%d/%d umi(total=%d allele=%d)", AlleleCount, Depth, UmiCounts.totalCount(), UmiCounts.alleleCount());
    }
}
