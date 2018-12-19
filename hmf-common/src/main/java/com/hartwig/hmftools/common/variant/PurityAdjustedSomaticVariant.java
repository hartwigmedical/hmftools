package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;

import org.jetbrains.annotations.NotNull;

public interface PurityAdjustedSomaticVariant extends SomaticVariant {

    double adjustedCopyNumber();

    double adjustedVAF();

    double minorAllelePloidy();

    default double ploidy() {
        return adjustedCopyNumber() * adjustedVAF();
    }

    default boolean biallelic() {
        return Doubles.lessOrEqual(adjustedCopyNumber(), 0) || Doubles.greaterOrEqual(ploidy(), adjustedCopyNumber() - 0.5);
    }

    @NotNull
    GermlineStatus germlineStatus();
}
