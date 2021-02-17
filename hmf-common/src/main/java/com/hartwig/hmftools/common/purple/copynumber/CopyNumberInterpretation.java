package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
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
        if (cnaDriver.driver() == DriverType.AMP) {
            return FULL_GAIN;
        }

        if (cnaDriver.driver() == DriverType.PARTIAL_AMP) {
            return PARTIAL_GAIN;
        }

        if (cnaDriver.driver() == DriverType.DEL) {
            return Doubles.greaterThan(cnaDriver.maxCopyNumber(), 0.5) ? PARTIAL_LOSS : FULL_LOSS;
        }

        throw new IllegalStateException("Driver not an AMP or DEL: " + cnaDriver);
    }
}