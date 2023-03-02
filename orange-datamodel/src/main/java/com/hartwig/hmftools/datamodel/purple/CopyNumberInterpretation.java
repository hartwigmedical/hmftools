package com.hartwig.hmftools.datamodel.purple;

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
}