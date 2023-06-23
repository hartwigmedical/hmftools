package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.String.format;

import com.hartwig.hmftools.ctdna.purity.variant.UmiTypeCounts;

public class GenotypeFragments
{
    public final String SampleName;
    public final int AlleleCount;
    public final int Depth;
    public final double QualTotal;
    public final UmiTypeCounts UmiCounts;

    public GenotypeFragments(final String sampleName, final int alleleCount, final int depth, final double qualTotal,
            final UmiTypeCounts umiCounts)
    {
        SampleName = sampleName;
        AlleleCount = alleleCount;
        Depth = depth;
        QualTotal = qualTotal;
        UmiCounts = umiCounts;
    }

    public double qualPerAlleleFragment()
    {
        return AlleleCount > 0 ? QualTotal / AlleleCount : 0;
    }

    public String toString() { return format("%d/%d umi(ref=%d allele=%d)",
            AlleleCount, Depth, UmiCounts.refTotal(), UmiCounts.alleleTotal()); }
}
