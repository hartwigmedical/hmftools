package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public interface StructuralVariantLegCopyNumber {

    byte orientation();

    @NotNull
    Optional<Double> leftCopyNumber();

    @NotNull
    Optional<Double> rightCopyNumber();

    default double adjustedCopyNumber() {
        if (orientation() == 1) {
            return leftCopyNumber().orElse(0D);
        } else {
            return rightCopyNumber().orElse(0D);
        }
    }

    default double adjustedCopyNumberChange() {
        double leftCopyNumber = leftCopyNumber().orElse(0D);
        double rightCopyNumber = rightCopyNumber().orElse(0D);

        return orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

}
