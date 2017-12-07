package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;

public interface PurityAdjustedSomaticVariantBuilder {

    PurityAdjustedSomaticVariant build();

    PurityAdjustedSomaticVariantBuilder adjustedVAF(double adjustedCopyNumber);

    PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(double adjustedVAF);

    PurityAdjustedSomaticVariantBuilder clonality(Clonality clonality);

    PurityAdjustedSomaticVariantBuilder lossOfHeterozygosity(boolean lossOfHeterozygosity);

    default PurityAdjustedSomaticVariantBuilder purityAdjustment(@NotNull PurityAdjuster purityAdjuster,
            @NotNull final PurpleCopyNumber copyNumber, @NotNull final AllelicDepth depth) {
        double adjustedCopyNumber = copyNumber.averageTumorCopyNumber();
        double adjustedVAF = purityAdjuster.purityAdjustedVAF(copyNumber.chromosome(), Math.max(0.001, adjustedCopyNumber), depth.alleleFrequency());
        double variantPloidy = adjustedCopyNumber * adjustedVAF;
        boolean loh = Doubles.lessOrEqual(adjustedCopyNumber, 0) || Doubles.greaterOrEqual(variantPloidy, adjustedCopyNumber - 0.5);
        return adjustedCopyNumber(adjustedCopyNumber).adjustedVAF(adjustedVAF).clonality(Clonality.UNKNOWN).lossOfHeterozygosity(loh);
    }

    default PurityAdjustedSomaticVariantBuilder clonality(@NotNull final PurityAdjuster purityAdjuster, @NotNull final PurpleCopyNumber copyNumber, @NotNull final AllelicDepth depth) {
        return clonality(Clonality.fromSample(purityAdjuster, copyNumber, depth));
    }
}
