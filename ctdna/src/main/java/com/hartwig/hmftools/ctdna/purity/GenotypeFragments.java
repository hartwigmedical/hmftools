package com.hartwig.hmftools.ctdna.purity;

public class GenotypeFragments
{
    public final String SampleName;
    public final int AlleleCount;
    public final int Depth;
    public final double QualTotal;

    public GenotypeFragments(final String sampleName, final int alleleCount, final int depth, final double qualTotal)
    {
        SampleName = sampleName;
        AlleleCount = alleleCount;
        Depth = depth;
        QualTotal = qualTotal;
    }

    public double qualPerAlleleFragment()
    {
        return AlleleCount > 0 ? QualTotal / AlleleCount : 0;
    }
}
