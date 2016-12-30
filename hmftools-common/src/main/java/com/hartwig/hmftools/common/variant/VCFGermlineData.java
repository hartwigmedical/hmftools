package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public class VCFGermlineData {

    @NotNull
    private final String refData;
    @NotNull
    private final String tumorData;
    @NotNull
    private final VCFType type;

    public VCFGermlineData(@NotNull final VCFType type, @NotNull final String refData,
            @NotNull final String tumorData) {
        this.refData = refData;
        this.tumorData = tumorData;
        this.type = type;
    }

    @NotNull
    public VCFType getType() {
        return type;
    }

    @NotNull
    public String getRefData() {
        return refData;
    }

    @NotNull
    public String getTumorData() {
        return tumorData;
    }
}
