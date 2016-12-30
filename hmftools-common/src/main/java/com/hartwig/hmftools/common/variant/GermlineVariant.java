package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public class GermlineVariant {

    @NotNull
    private final VariantType type;
    @NotNull
    private final String refData;
    @NotNull
    private final String tumorData;

    public GermlineVariant(@NotNull final VariantType type, @NotNull final String refData,
            @NotNull final String tumorData) {
        this.type = type;
        this.refData = refData;
        this.tumorData = tumorData;
    }

    @NotNull
    public VariantType type() {
        return type;
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
