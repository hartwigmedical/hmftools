package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.String.format;

import com.hartwig.hmftools.common.variant.AllelicDepth;

public class WeightedPloidy extends AllelicDepth
{
    public double Ploidy;
    public double Weight;

    public WeightedPloidy(final int totalReadCount, final int alleleReadCount, final double ploidy, final double weight)
    {
        super(totalReadCount, alleleReadCount);
        Ploidy = ploidy;
        Weight = weight;
    }

    public double ploidy() { return Ploidy; }
    public double weight() { return Weight; }

    public String toString() { return format("frags(%s) ploidy(%.2f) weight(%.2f)", super.toString(), Ploidy, Weight); }

    public static WeightedPloidy create(double ploidy, int alleleReadCount, int totalReadCount)
    {
        return new WeightedPloidy(totalReadCount, alleleReadCount, ploidy, 1);
    }
}
