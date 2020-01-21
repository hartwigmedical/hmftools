package com.hartwig.hmftools.common.variant.tml;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum TumorMutationalStatus {
    HIGH,
    LOW,
    UNKNOWN;

    private static final double TMB_CUTOFF = 10;
    private static final double TML_CUTOFF = 130;

    @NotNull
    public static TumorMutationalStatus fromBurdenPerMb(double burdenPerMb) {
        return Doubles.greaterThan(burdenPerMb, TMB_CUTOFF) ? HIGH : LOW;
    }

    @NotNull
    public static TumorMutationalStatus fromLoad(double load) {
        return Doubles.greaterThan(load, TML_CUTOFF) ? HIGH : LOW;
    }

}
