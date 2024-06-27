package com.hartwig.hmftools.purple.fitting;

import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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

    public static WeightedPloidy create(double ploidy, int alleleReadCount, int totalReadCount)
    {
        return new WeightedPloidy(totalReadCount, alleleReadCount, ploidy, 1);
    }
}
