package com.hartwig.hmftools.common.variant;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class SomaticVariant {

    @NotNull
    private final VariantType type;
    @NotNull
    private final List<String> callers;
    private final double alleleFrequency;
    private final boolean isDBSNP;
    private final boolean isCOSMIC;

    SomaticVariant(@NotNull final VariantType type, @NotNull final List<String> callers, final double alleleFrequency,
            final boolean isDBSNP, final boolean isCOSMIC) {
        this.type = type;
        this.callers = callers;
        this.alleleFrequency = alleleFrequency;
        this.isDBSNP = isDBSNP;
        this.isCOSMIC = isCOSMIC;
    }

    @NotNull
    public VariantType type() {
        return type;
    }

    @NotNull
    public List<String> callers() {
        return callers;
    }

    public long callerCount() {
        return callers.size();
    }

    public double alleleFrequency() {
        return alleleFrequency;
    }

    public boolean isDBSNP() {
        return isDBSNP;
    }

    public boolean isCOSMIC() {
        return isCOSMIC;
    }
}
