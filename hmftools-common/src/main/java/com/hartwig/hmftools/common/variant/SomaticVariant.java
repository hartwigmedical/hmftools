package com.hartwig.hmftools.common.variant;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class SomaticVariant implements Variant {

    @NotNull
    private final VariantType type;
    @NotNull
    private final String filter;
    @NotNull
    private final String chromosome;
    private final int position;

    @NotNull
    private final List<String> callers;
    private final double alleleFrequency;
    private final boolean isDBSNP;
    private final boolean isCOSMIC;

    SomaticVariant(@NotNull final VariantType type, @NotNull final String filter, @NotNull final String chromosome,
            final int position, @NotNull final List<String> callers, final double alleleFrequency,
            final boolean isDBSNP, final boolean isCOSMIC) {
        this.type = type;
        this.filter = filter;
        this.chromosome = chromosome;
        this.position = position;
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
    public String filter() {
        return filter;
    }

    @NotNull
    String chromosome() {
        return chromosome;
    }

    int position() {
        return position;
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
