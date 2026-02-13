package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SageVcfTags;

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

    private final List<FilterReason> mFilterReasons;
    private final List<FilterReason> mDualFilterReasons;

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
        mFilterReasons = Lists.newArrayList();
        mDualFilterReasons = Lists.newArrayList();
    }

    public double qualPerAlleleFragment() { return AlleleCount > 0 ? QualTotal / (double)AlleleCount : 0; }

    public double vaf() { return Depth > 0 ? AlleleCount / (double)Depth : 0; }

    public void setBqrErrorRate(double errorRate) { mBqrErrorRate = errorRate; }
    public double bqrErrorRate() { return mBqrErrorRate; }

    public void markOutlier()
    {
        mFilterReasons.add(FilterReason.OUTLIER);
    }

    public double averageReadDistance()
    {
        Object value = GenotypeData.getExtendedAttribute(SageVcfTags.AVG_EDGE_DISTANCE_PERC, null);

        if(value == null)
            return -1;

        String[] aedCounts = value.toString().split(LIST_SEPARATOR, 2);
        return aedCounts.length == 2 ? Double.parseDouble(aedCounts[1]) : -1;
    }

    public boolean isFiltered() { return !mFilterReasons.isEmpty(); }
    public void addFilterReason(final FilterReason filterReason) { mFilterReasons.add(filterReason); }
    public List<FilterReason> filterReasons() { return mFilterReasons; }

    public boolean isDualFiltered() { return !mDualFilterReasons.isEmpty(); }
    public void addDualFilterReason(final FilterReason filterReason) { mDualFilterReasons.add(filterReason); }
    public List<FilterReason> dualFilterReasons() { return mDualFilterReasons; }

    public String toString()
    {
        return format("%d/%d umi(total=%d allele=%d)", AlleleCount, Depth, UmiCounts.totalCount(), UmiCounts.alleleCount());
    }
}
