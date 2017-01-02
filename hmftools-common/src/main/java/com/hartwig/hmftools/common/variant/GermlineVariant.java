package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public class GermlineVariant implements Variant {

    @NotNull
    private final VariantType type;
    @NotNull
    private final String filter;
    @NotNull
    private final String refData;
    @NotNull
    private final String tumorData;

    public GermlineVariant(@NotNull final VariantType type, @NotNull final String filter,
            @NotNull final String refData, @NotNull final String tumorData) {
        this.type = type;
        this.filter = filter;
        this.refData = refData;
        this.tumorData = tumorData;
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
    public String refData() {
        return refData;
    }

    @NotNull
    public String tumorData() {
        return tumorData;
    }
}
