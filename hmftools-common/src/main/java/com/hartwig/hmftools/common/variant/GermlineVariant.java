package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GermlineVariant implements Variant {

    @NotNull
    private final VariantType type;
    @NotNull
    private final String filter;
    @NotNull
    private final GermlineSampleData refData;
    @Nullable
    private final GermlineSampleData tumorData;

    public GermlineVariant(@NotNull final VariantType type, @NotNull final String filter,
            @NotNull final GermlineSampleData refData, @Nullable final GermlineSampleData tumorData) {
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
    public GermlineSampleData refData() {
        return refData;
    }

    @Nullable
    public GermlineSampleData tumorData() {
        return tumorData;
    }
}
