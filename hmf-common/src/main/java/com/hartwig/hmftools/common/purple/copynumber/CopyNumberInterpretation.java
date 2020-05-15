package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberInterpretation {
    GAIN("gain"),
    FULL_LOSS("full loss"),
    PARTIAL_LOSS("partial loss");

    @NotNull
    private final String text;

    CopyNumberInterpretation(@NotNull final String text) {
        this.text = text;
    }

    @NotNull
    public String text() {
        return text;
    }

    @NotNull
    public static CopyNumberInterpretation fromCNADriver(@NotNull DriverCatalog cnaDriver) {
        if (cnaDriver.driver() == DriverType.AMP) {
            return GAIN;
        }

        if (cnaDriver.driver() == DriverType.DEL) {
            return Doubles.greaterThan(cnaDriver.maxCopyNumber(), 0.5) ? PARTIAL_LOSS : FULL_LOSS;
        }

        throw new IllegalStateException("Driver not an AMP or DEL: " + cnaDriver);
    }
}