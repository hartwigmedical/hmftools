package com.hartwig.hmftools.purple.fitting;

import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface WeightedPloidy extends AllelicDepth
{
    double ploidy();

    double weight();

    @NotNull
    static WeightedPloidy create(double ploidy, int alleleReadCount, int totalReadCount)
    {
        return ModifiableWeightedPloidy.create()
                .setPloidy(ploidy)
                .setWeight(1)
                .setAlleleReadCount(alleleReadCount)
                .setTotalReadCount(totalReadCount);
    }
}
