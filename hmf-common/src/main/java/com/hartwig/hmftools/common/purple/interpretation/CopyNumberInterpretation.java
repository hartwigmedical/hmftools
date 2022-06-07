package com.hartwig.hmftools.common.purple.interpretation;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberInterpretation {
    FULL_GAIN("full gain"),
    PARTIAL_GAIN("partial gain"),
    FULL_LOSS("full loss"),
    PARTIAL_LOSS("partial loss");

    @NotNull
    private final String display;

    CopyNumberInterpretation(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

    @NotNull
    public static CopyNumberInterpretation fromCNADriver(@NotNull DriverCatalog cnaDriver) {
        switch (cnaDriver.driver()) {
            case AMP:
                return FULL_GAIN;
            case PARTIAL_AMP:
                return PARTIAL_GAIN;
            case DEL:
                return Doubles.greaterThan(cnaDriver.maxCopyNumber(), 0.5) ? PARTIAL_LOSS : FULL_LOSS;
            default:
                throw new IllegalStateException("Driver not an AMP or DEL: " + cnaDriver);
        }
    }
}