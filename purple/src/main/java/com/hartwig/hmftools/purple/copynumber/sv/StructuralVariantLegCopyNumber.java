package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;

public interface StructuralVariantLegCopyNumber extends StructuralVariantLeg {

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
}
