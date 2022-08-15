package com.hartwig.hmftools.common.purple.loader;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberInterpretation {
    FULL_GAIN,
    PARTIAL_GAIN,
    FULL_LOSS,
    PARTIAL_LOSS;

    @NotNull
    public String display() {
        return this.toString().toLowerCase().replaceAll("_", " ");
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