package com.hartwig.hmftools.common.variant;

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
        double adjustedVAF = purityAdjuster.purityAdjustedVAF(Math.max(0.001, adjustedCopyNumber), depth.alleleFrequency());
        return adjustedCopyNumber(adjustedCopyNumber).adjustedVAF(adjustedVAF).clonality(Clonality.UNKNOWN).lossOfHeterozygosity(false);
    }
}
